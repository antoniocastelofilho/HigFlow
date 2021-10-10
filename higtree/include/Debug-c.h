#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <assert.h>
#include <sys/timeb.h>
#ifdef DEBUG

static struct timeb __lasttime;
static struct timeb __thistime;
#define __DEBUG_STACK 100
static int __debug_flag[__DEBUG_STACK] = {1};
static int __debug_top = 0;

#define debugfd stdout
//! Prints the time when this line was executed, if DEBUG is defined
#define DEBUG_TIME {ftime(&__thistime); __lasttime = __thistime; if (__debug_flag[__debug_top]) {fprintf(debugfd, "%s:%d: time %ld %d\n", __FILE__, __LINE__, __thistime.time, __thistime.millitm);}}

//! Prints how many milliseconds ellapsed between two executions of this macro, if DEBUG is defined
#define DEBUG_DIFF_TIME {__lasttime = __thistime; ftime(&__thistime); if (__debug_flag[__debug_top]) {fprintf(debugfd, "%s:%d: elapsed time %ld miliseconds\n", __FILE__, __LINE__, (1000 * __thistime.time + __thistime.millitm) - (1000 * __lasttime.time + __lasttime.millitm));}}

#define DEBUG_DIFF_TIME_START {ftime(&__lasttime);}

#define DEBUG_DIFF_TIME_FINISH(label) {ftime(&__thistime); fprintf(debugfd, "%s: %.03f seconds\n", label, (__thistime.time + __thistime.millitm / 1000.0) - (__lasttime.time + __lasttime.millitm / 1000.0));}

//! Executes the command x, if DEBUG is defined
#define DEBUG_EXEC(x) do { if (__debug_flag[__debug_top]) {DEBUG_PASS; x; }} while(0)

//! Prints the filename and line of the command, if DEBUG is defined
#define DEBUG_PASS {if (__debug_flag[__debug_top]) {fprintf(debugfd, "%s:%d\n", __FILE__, __LINE__);}}

//! Prints the message x, if DEBUG is defined
#define DEBUG_WARNING(x) {fprintf(debugfd, x);}

//! Prints the contents of variable x, using printf format code f, if DEBUG is defined
#define DEBUG_INSPECT(x, f) {if (__debug_flag[__debug_top]) {fprintf(debugfd, "%s:%d %s = " #f "\n", __FILE__, __LINE__, #x, x);}}

#define DEBUG_PUSH(b) {if (__debug_top < __DEBUG_STACK) { __debug_flag[++__debug_top] = (b);} }

#define DEBUG_POP {if (__debug_top > 0) {__debug_top--;} }

#define DEBUG_ASSERT(x) {if (!(x)) {DEBUG_PASS; fprintf(debugfd, "assertion failed: %s\n", #x); *(int *)NULL = 0;}}

#else

#define DEBUG_TIME
#define DEBUG_DIFF_TIME ;
#define DEBUG_EXEC(x)
#define DEBUG_PASS
#define DEBUG_WARNING(x)
#define DEBUG_INSPECT(x, f) ;
#define DEBUG_PUSH(b) ;
#define DEBUG_POP ;
#define DEBUG_DIFF_TIME_START
#define DEBUG_DIFF_TIME_FINISH(x)
#define DEBUG_ASSERT(x) ;

#endif

//! Checks if retcode is an MPI error, and prints a message in debug output
//! Ensures retcode is evaluated once, and only once.
#define CHECK_MPI_ERR(retcode) do{ \
	int rcode = (retcode); \
	if(rcode != MPI_SUCCESS) { \
		char msg[MPI_MAX_ERROR_STRING]; \
		int len; \
		MPI_Error_string(rcode, msg, &len); \
		printf("MPI ERROR on %s:%d, (%s):\n   %s \n", \
			__FILE__, __LINE__, __func__, msg); \
	} \
} while(0)

#endif
