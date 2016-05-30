#ifndef SSTRING_INCLUDED
#define SSTRING_INCLUDED

#include <cstdio>
#include <string>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>

// RKB -- 5/29/08
#include <iterator>

#include "util/misc.h"

// default maximum length for sstring::getline
#define DART_MAX_LINE_SIZE 12345678

namespace fsa {

  // base class for sstring
  class string_with_streambuf : public std::basic_string<char>
    {
    protected:
      // streambuf prototype
      class DART_string_streambuf : public std::streambuf
	{
	protected:
	  int overflow(int c);
	public:
	  string_with_streambuf& owner;
	  DART_string_streambuf (string_with_streambuf& s) : owner(s) { }
	};
  
      // streambuf member
      DART_string_streambuf sbuf;

    public:
      // constructors
    string_with_streambuf() : std::basic_string<char>(), sbuf(*this) { }
      string_with_streambuf (long n, const char value) : std::basic_string<char>(n,value), sbuf(*this) { }
      string_with_streambuf (const std::basic_string<char>& x) : std::basic_string<char>(x), sbuf(*this) { }
      string_with_streambuf (const string_with_streambuf& x) : std::basic_string<char>(x), sbuf(*this) { }
#ifdef __STL_MEMBER_TEMPLATES
      template <class InputIterator>
	string_with_streambuf(InputIterator first, InputIterator last) : std::basic_string<char>(first, last), sbuf(*this) { }
#else /* __STL_MEMBER_TEMPLATES */
    string_with_streambuf(const_iterator first, const_iterator last) : std::basic_string<char>(first, last), sbuf(*this) { }
#endif /* __STL_MEMBER_TEMPLATES */
      string_with_streambuf (const char* x) : std::basic_string<char>(x), sbuf(*this) { }
      string_with_streambuf (const char* x, size_type n) : std::basic_string<char>(x,n), sbuf(*this) { }
    };

  // sstring: DART string class
  class sstring : public string_with_streambuf, public std::ostream
  {
  public:
  sstring() : string_with_streambuf(), std::ostream(&sbuf) {}
    sstring (size_type n, const char value = '\0') : string_with_streambuf(n,value), std::ostream(&sbuf) {}
    sstring (int n, const char value = '\0') : string_with_streambuf(n,value), std::ostream(&sbuf) {}
    sstring (long n, const char value = '\0') : string_with_streambuf(n,value), std::ostream(&sbuf) {}
    sstring (const std::basic_string<char>& x) : string_with_streambuf(x), std::ostream(&sbuf) {}
    sstring (const sstring& x) : string_with_streambuf(x), std::ostream(&sbuf) {}

#ifdef __STL_MEMBER_TEMPLATES
    template <class InputIterator>
      sstring(InputIterator first, InputIterator last) : string_with_streambuf(first, last), std::ostream(&sbuf) {}
#else /* __STL_MEMBER_TEMPLATES */
  sstring(const_iterator first, const_iterator last) : string_with_streambuf(first, last), std::ostream(&sbuf) {}
#endif /* __STL_MEMBER_TEMPLATES */

    sstring (const char* x) : string_with_streambuf(x), std::ostream(&sbuf) {}
    sstring (const char* x, size_type n) : string_with_streambuf(x,n), std::ostream(&sbuf) {}

    void read_entire_file (const char* filename);

    iterator end() { return ((std::basic_string<char>*) this) -> end(); }
    const_iterator end() const { return ((const std::basic_string<char>*) this) -> end(); }

    bool operator==(const std::basic_string<char>& x) const { return ((const std::basic_string<char>&) *this) == x; }
    bool operator==(const char* x) const { return ((const std::basic_string<char>&) *this) == std::basic_string<char> (x); }

    bool operator<(const std::basic_string<char>& x) const { return ((const std::basic_string<char>&) *this) < x; }

    sstring& operator=(const sstring& x) { ((std::basic_string<char>&) *this) = (const std::basic_string<char>&) x; return *this; }
    sstring& operator=(const std::basic_string<char>& x) { ((std::basic_string<char>&) *this) = x; return *this; }
    sstring& operator=(const char* x) { ((std::basic_string<char>&) *this) = x; return *this; }

    void push_back (char c) { append (1, c); }
    void clear() { erase(); }
    void to_lower();
    void to_upper();

    std::vector<sstring> split (const char* split_chars = " \t\n", bool skip_empty_fields = 1, int max_fields = 0) const;
    static sstring join (const std::vector<sstring>& v, const char* sep = " ");

    friend std::ostream& operator<<(std::ostream&o, const sstring& s) { return o << s.c_str(); }

    sstring& getline (std::istream& is, size_type max_size = DART_MAX_LINE_SIZE);
  
    char front() { char c = '\0'; if (size() > 0) c = *begin(); return c; }
    char back() { char c = '\0'; if (size() > 0) c = *(end()-1); return c; }
    char chop() { char c = '\0'; if (size() > 0) { c = back(); erase (end() - 1); } return c; }
    char chomp (const char chomp_char = '\n') { return back() == chomp_char ? chop() : '\0'; }

    int to_int() const { return atoi (c_str()); }
    double to_double() const { return atof (c_str()); }

    // strict versions of to_int() and to_double(), that do checking
    int to_int_strict (const char* err_prefix = "Not an integer constant") const;
    int to_nonneg_int_strict (const char* err_prefix = "Not a non-negative integer constant") const;
    double to_double_strict (const char* err_prefix = "Not a numeric constant") const;
    double to_nonneg_double_strict (const char* err_prefix = "Not a non-negative numeric constant") const;

    sstring substr (int start, int len) const;

  };

  // char_string converts non-printable chars into octal
  class char_string : public sstring
  {
  private:
    void add_char (char c);
  public:
    char_string (char c);     
    char_string (const char* s);
  };

}

#endif
