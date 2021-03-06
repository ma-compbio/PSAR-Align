/***************************Copyright-DO-NOT-REMOVE-THIS-LINE**
  *
  * Condor Software Copyright Notice
  * Copyright (C) 1990-2004, Condor Team, Computer Sciences Department,
  * University of Wisconsin-Madison, WI.
  *
  * This source code is covered by the Condor Public License, which can
  * be found in the accompanying LICENSE.TXT file, or online at
  * www.condorproject.org.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  * AND THE UNIVERSITY OF WISCONSIN-MADISON "AS IS" AND ANY EXPRESS OR
  * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  * WARRANTIES OF MERCHANTABILITY, OF SATISFACTORY QUALITY, AND FITNESS
  * FOR A PARTICULAR PURPOSE OR USE ARE DISCLAIMED. THE COPYRIGHT
  * HOLDERS AND CONTRIBUTORS AND THE UNIVERSITY OF WISCONSIN-MADISON
  * MAKE NO MAKE NO REPRESENTATION THAT THE SOFTWARE, MODIFICATIONS,
  * ENHANCEMENTS OR DERIVATIVE WORKS THEREOF, WILL NOT INFRINGE ANY
  * PATENT, COPYRIGHT, TRADEMARK, TRADE SECRET OR OTHER PROPRIETARY
  * RIGHT.
  *
  ****************************Copyright-DO-NOT-REMOVE-THIS-LINE**/
/*
Chirp C Client

This public domain software is provided "as is".  See the Chirp License
for details.
*/

#include "chirp_protocol.h"
#include "chirp_client.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include <unistd.h>
#include <sys/errno.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <netinet/tcp.h>

static int tcp_connect( const char *host, int port );
static void chirp_fatal_request( const char *name );
static void chirp_fatal_response();
static int get_result( FILE *s );
static int convert_result( int response );
static int simple_command(struct chirp_client *c,char const *fmt,...);
static void vsprintf_chirp(char *command,char const *fmt,va_list args);
static char const *read_url_param(char const *url,char *buffer,size_t length);


struct chirp_client {
	FILE *stream;
};

/*
  chirp_client_connect_url()

  Sets path_part to the position of the start of the path information
  in the URL and connects to the specified Chirp server.

  URL format:
    chirp:host.name:port/path   (absolute path)
    chirp:host.name:port./path  (relative path)
    chirp:/path                 (absolute path)
    chirp:./path                (relative path)
    chirp:path                  (sloppy relative path)

  Note that the initial part of a sloppy relative path can be confused
  for a host:port specification if it happens to look like one.  Example:
  'chirp:strange:1/file'.  For this reason, it is better to use one of
  the non-sloppy formats.

  In all of the above URL formats, any number of extra connection
  parameters may be inserted before the path part (including  the leading
  '/' or '.').  These are of the form ";parameter=value" where parameter
  and value are encoded using the standard Mime type
  application/x-www-form-url encoded, just like the parameters in a typical
  HTTP GET CGI request, but using ';' instead of '&' as the delimiter.

  At this time, no connection parameters are defined, and any that
  are supplied are simply ignored.

*/

struct chirp_client *
chirp_client_connect_url( const char *url, const char **path_part)
{
	struct chirp_client *client;
	char const *str;
	char *host = NULL;
	int port = 0;

	if(strncmp(url,"chirp:",6)) {
		//bare file name
		*path_part = url;
		return chirp_client_connect_default();
	}

	url += 6; // "chirp:"

	if(*url != '/' && *url != '\\' && *url != ';' && *url != '.' \
	   && (str = strchr(url,':')))
	{
		char *end;
		port = strtol(str+1,&end,10);
		if(port && end > str+1 && 
		   (*end == '\0' || *end == '/' || *end == '\\' ||
		    *end == '.' || *end == ';'))
		{
			//URL form chirp:host.name:port...
			//Note that we try to avoid getting here on a "sloppy"
			//relative path that happens to contain a ':' but
			//which is not followed by a valid port/path.

			host = (char *)malloc(str-url+1);
			strncpy(host,url,str-url);
			host[str-url] = '\0';

			url = end;
		}
	}

	while(*url == ';') { //parse connection parameters
		char param[CHIRP_LINE_MAX];
		char value[CHIRP_LINE_MAX];

		url = read_url_param(++url,param,sizeof(param));
		if(!url) {
			errno = EINVAL;
			return NULL;
		}

		if(*url == '=') {
			url = read_url_param(++url,value,sizeof(value));
			if(!url) {
				errno = EINVAL;
				return NULL;
			}
		}
		else *value = '\0';

		//No connection parameters are defined at this time!
		//Handle them here when they are defined.
	}

	*path_part = url;

	if(!host) { //URL must be in form 'chirp:path'
		client = chirp_client_connect_default();
	}
	else {
		client = chirp_client_connect(host,port);
	}

	free(host);
	return client;
}

struct chirp_client *
chirp_client_connect_default()
{
	FILE *file;
	int fields;
	int save_errno;
	struct chirp_client *client;
	char host[CHIRP_LINE_MAX];
	char cookie[CHIRP_LINE_MAX];
	int port;
	int result;

	file = fopen("chirp.config","r");
	if(!file) return 0;

	fields = fscanf(file,"%s %d %s",host,&port,cookie);
	fclose(file);

	if(fields!=3) {
		errno = EINVAL;
		return 0;
	}

	client = chirp_client_connect(host,port);
	if(!client) return 0;

	result = chirp_client_cookie(client,cookie);
	if(result!=0) {
		save_errno = errno;
		chirp_client_disconnect(client);
		errno = save_errno;
		return 0;
	}

	return client;
}

struct chirp_client *
chirp_client_connect( const char *host, int port )
{
	struct chirp_client *c;
	int save_errno;
	int fd;

	c = (chirp_client *) malloc(sizeof(*c));
	if(!c) return 0;

	fd = tcp_connect(host,port);
	if(fd<0) {
		save_errno = errno;
		free(c);
		errno = save_errno;
		return 0;
	}

	c->stream = fdopen(fd,"w+");
	if(!c->stream) {
		save_errno = errno;
		close(fd);
		free(c);
		errno = save_errno;
		return 0;
	}

	return c;
}

void
chirp_client_disconnect( struct chirp_client *c )
{
	fclose(c->stream);
	free(c);
}

int
chirp_client_cookie( struct chirp_client *c, const char *cookie )
{
	return simple_command(c,"cookie %s\n",cookie);
}

int
chirp_client_login( struct chirp_client *c, const char *name, const char *password )
{
	return simple_command(c,"login %s %s\n",name,password);
}

int
chirp_client_lookup( struct chirp_client *c, const char *logical_name, char **url )
{
	int result;
	int actual;

	result = simple_command(c,"lookup %s\n",logical_name);
	if(result>0) {
		*url = (char *)malloc(result);
		if(*url) {
			actual = fread(*url,1,result,c->stream);
			if(actual!=result) chirp_fatal_request("lookup");
		} else {
			chirp_fatal_request("lookup");
		}
	}

	return result;
}

int
chirp_client_constrain( struct chirp_client *c, const char *expr)
{
	return simple_command(c,"constrain %s\n",expr);
}

int
chirp_client_get_job_attr( struct chirp_client *c, const char *name, char **expr )
{
	int result;
	int actual;

	result = simple_command(c,"get_job_attr %s\n",name);
	if(result>0) {
		*expr = (char *)malloc(result);
		if(*expr) {
			actual = fread(*expr,1,result,c->stream);
			if(actual!=result) chirp_fatal_request("get_job_attr");
		} else {
			chirp_fatal_request("get_job_attr");
		}
	}

	return result;
}

int
chirp_client_set_job_attr( struct chirp_client *c, const char *name, const char *expr )
{
	return simple_command(c,"set_job_attr %s %s\n",name,expr);
}

int
chirp_client_open( struct chirp_client *c, const char *path, const char *flags, int mode )
{
	return simple_command(c,"open %s %s %d\n",path,flags,mode);
}

int
chirp_client_close( struct chirp_client *c, int fd )
{
	return simple_command(c,"close %d\n",fd);
}

int
chirp_client_read( struct chirp_client *c, int fd, void *buffer, int length )
{
	int result;
	int actual;

	result = simple_command(c,"read %d %d\n",fd,length);

	if( result>0 ) {
		actual = fread(buffer,1,result,c->stream);
		if(actual!=result) chirp_fatal_request("read");
	}

	return result;
}

int
chirp_client_write( struct chirp_client *c, int fd, const void *buffer, int length )
{
	int actual;
	int result;

	result = fprintf(c->stream,"write %d %d\n",fd,length);
	if(result<0) chirp_fatal_request("write");

	result = fflush(c->stream);
	if(result<0) chirp_fatal_request("write");

	actual = fwrite(buffer,1,length,c->stream);
	if(actual!=length) chirp_fatal_request("write");

	return convert_result(get_result(c->stream));
}

int
chirp_client_unlink( struct chirp_client *c, const char *path )
{
	return simple_command(c,"unlink %s\n",path);
}

int
chirp_client_rename( struct chirp_client *c, const char *oldpath, const char *newpath )
{
	return simple_command(c,"rename %s %s\n",oldpath,newpath);
}
int
chirp_client_fsync( struct chirp_client *c, int fd )
{
	return simple_command(c,"fsync %d\n",fd);
}

int
chirp_client_lseek( struct chirp_client *c, int fd, int offset, int whence )
{
	return simple_command(c,"lseek %d %d %d\n",fd,offset,whence);
}

int
chirp_client_mkdir( struct chirp_client *c, char const *name, int mode )
{
	return simple_command(c,"mkdir %s %d\n",name,mode);
}

int
chirp_client_rmdir( struct chirp_client *c, char const *name )
{
	return simple_command(c,"rmdir %s\n",name);
}


static int
convert_result( int result )
{
	if(result>=0) {
		return result;
	} else {
		switch(result) {
			case CHIRP_ERROR_NOT_AUTHENTICATED:
			case CHIRP_ERROR_NOT_AUTHORIZED:
				errno = EACCES;
				break;
			case CHIRP_ERROR_DOESNT_EXIST:
				errno = ENOENT;
				break;
			case CHIRP_ERROR_ALREADY_EXISTS:
				errno = EEXIST;
				break;
			case CHIRP_ERROR_TOO_BIG:
				errno = EFBIG;
				break;
			case CHIRP_ERROR_NO_SPACE:
				errno = ENOSPC;
				break;
			case CHIRP_ERROR_NO_MEMORY:
				errno = ENOMEM;
				break;
			case CHIRP_ERROR_INVALID_REQUEST:
				errno = EINVAL;
				break;
			case CHIRP_ERROR_TOO_MANY_OPEN:
				errno = EMFILE;
				break;
			case CHIRP_ERROR_BUSY:
				errno = EBUSY;
				break;
			case CHIRP_ERROR_TRY_AGAIN:
				errno = EAGAIN;
				break;
			case CHIRP_ERROR_UNKNOWN:
				chirp_fatal_response();
				break;
		}
		return -1;
	}
}

static int
get_result( FILE *s )
{
	char line[CHIRP_LINE_MAX];
	char *c;
	int result;
	int fields;

	c = fgets(line,CHIRP_LINE_MAX,s);
	if(!c) chirp_fatal_response();

	fields = sscanf(line,"%d",&result);
	if(fields!=1) chirp_fatal_response();

#ifdef CHIRP_DEBUG
	fprintf(stderr,"chirp received: %s\n",line);
#endif

	return result;
}

static void
chirp_fatal_request( const char *name )
{
	fprintf(stderr,"chirp: couldn't %s: %s\n",name,strerror(errno));
	abort();
}

static
void chirp_fatal_response()
{
	fprintf(stderr,"chirp: couldn't get response from server: %s\n",strerror(errno));
	abort();
}

static int
tcp_connect( const char *host, int port )
{
	struct hostent *h;
	struct sockaddr_in address;
	int success;
	int fd;

	h = gethostbyname(host);
	if(!h) return -1;

	address.sin_port = htons(port);
	address.sin_family = h->h_addrtype;
	memcpy(&address.sin_addr.s_addr,h->h_addr_list[0],sizeof(address.sin_addr.s_addr));

	fd = socket( AF_INET, SOCK_STREAM, 0 );
	if(fd<0) return -1;

	success = connect( fd, (struct sockaddr *) &address, sizeof(address) );
	if(success<0) {
		close(fd);
		return -1;
	}

	return fd;
}

/*
  vsprintf_chirp -- simple sprintf capabilities with character escaping

  The following format characters are interpreted:

  %d -- decimal
  %s -- word (whitespace is escaped)
  %% -- output %
*/

void
vsprintf_chirp(char *command,char const *fmt,va_list args)
{
	char       *c;
	char const *f;

	c = command;
	f = fmt;
	while(*f) {
		if(*f == '%') {
			switch(*(++f)) {
			case 'd':
				f++;
				sprintf(c,"%d",va_arg(args,int));
				c += strlen(c);
				break;
			case 's': {
				char const *w = va_arg(args,char const *);
				f++;
				while(*w) {
					switch(*w) {
					case ' ':
					case '\t':
					case '\n':
					case '\r':
					case '\\':
						*(c++) = '\\';
						/*fall through*/
					default:
						*(c++) = *(w++);
					}
				}
				break;
			}
			case '%':
				*(c++) = *(f++);
				break;
			default:
				chirp_fatal_request(f);
			}
		} else {
			*(c++) = *(f++);
		}
	}
	*(c++) = '\0';
}

int
simple_command(struct chirp_client *c,char const *fmt,...)
{
	int     result;
	char    command[CHIRP_LINE_MAX];
	va_list args;

	va_start(args,fmt);
	vsprintf_chirp(command,fmt,args);
	va_end(args);

#ifdef DEBUG_CHIRP
	fprintf(stderr,"chirp sending: %s",command);
#endif

	result = fputs(command,c->stream);



	if(result < 0) chirp_fatal_request(fmt);

	result = fflush(c->stream);
	if(result < 0) chirp_fatal_request(fmt);

	return convert_result(get_result(c->stream));
}

char const *
read_url_param(char const *url,char *buffer,size_t length)
{
	size_t bufpos = 0;

	while(*url != '\0' && *url != '.' && *url != '/'
	      && *url != '=' && *url != ';' && *url != '\\')
	{
		if(bufpos >= length) return NULL;	

		switch(*url) {
		case '+':
			buffer[bufpos++] = ' ';
			break;
		case '%': { //form-url-encoded escape sequence
			//following two characters are hex digits
			char d = tolower(*(++url));

			if(d >= '0' && d <= '9') d -= '0';
			else if(d >= 'a' && d <= 'f') d = d - 'a' + 0xA;
			else return NULL; //invalid hex digit

			buffer[bufpos] = d<<4;

			d = tolower(*(++url));

			if(d >= '0' && d <= '9') d -= '0';
			else if(d >= 'a' && d <= 'f') d = d - 'a' + 0xA;
			else return NULL; //invalid hex digit

			buffer[bufpos++] |= d;

			break;
		}
		default:
			buffer[bufpos++] = *url;
			break;
		}

		url++;
	}

	if(bufpos >= length) return NULL;
	buffer[bufpos] = '\0';

	return url;
}
