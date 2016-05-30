
// This class implements all functionality outside of ANSI C/C++
// All the ifdefery should be located here, and every class
// above this should either call into MWSystem, or use
// portable, ANSI code.

class MWSystem {
 public:
	static int getpid();
	static void gethostname(char *, int);
	static double gettimeofday();
	static double getcputime();
	static void sleep(int);
	static void mkdir(char *);
};
