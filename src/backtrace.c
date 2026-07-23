//
// This file implements a symbolic backtrace when a program aborts
// Especially when we lanuch bigstick on hopper or edison we may
// not want to do another run to find the failing routine.
// Ken McElvain  23Oct2013
// 
#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>
#include <signal.h>
#include <sys/resource.h>

//
// This part is a real hack.  Apparently Apple decided that no one needs
// to know where a signal handler to going to return to or any other
// register values.
// However, the sigcontext pointer is still there.  The problen is to know
// how it is layed out.  I managed to basically find them in /usr/include/sys/_structs.h
// but the 64 bit version is not complete and leaves the last field off of ucontext.
// However, the data is still there.   I hope Apple makes this work or fixes backtrace.
//
// Something similar should work for linux
//

#ifdef __APPLE__
//#define	_STRUCT_X86_THREAD_STATE64	struct __darwin_x86_thread_state64
struct thread_state64
{
	__uint64_t	__rax;
	__uint64_t	__rbx;
	__uint64_t	__rcx;
	__uint64_t	__rdx;
	__uint64_t	__rdi;
	__uint64_t	__rsi;
	__uint64_t	__rbp;
	__uint64_t	__rsp;
	__uint64_t	__r8;
	__uint64_t	__r9;
	__uint64_t	__r10;
	__uint64_t	__r11;
	__uint64_t	__r12;
	__uint64_t	__r13;
	__uint64_t	__r14;
	__uint64_t	__r15;
	__uint64_t	__rip;
	__uint64_t	__rflags;
	__uint64_t	__cs;
	__uint64_t	__fs;
	__uint64_t	__gs;
};

typedef struct my_ucontext {
	unsigned long pad[6];
	struct thread_state64 regs;
} my_ucontext_t;

typedef struct my_sig_ucontext {
	__uint64_t  pad; 
	my_ucontext_t   *uc_link;
} my_sig_ucontext_t;
#else
typedef void *my_sig_ucontext_t;
#endif

/* we want to avoid malloc happening in a signal handler, just allocate static space */
#define MAXFRAMES 64
static void *addrlist[MAXFRAMES+1];

static void myfputs(char *str, FILE *ofp) {
	fputs(str, ofp);
	fflush(ofp);
}

struct stack_frame {
	struct stack_frame *next;
	void *ret;
};

#if 0
//
// Do not delete this code.  It is an example of how to extract a stack frame on x86_64
// http://www.opensource.apple.com/source/openmpi/openmpi-5/openmpi/opal/mca/backtrace/ 
//
static int  mybacktrace(void **abuf, int maxframes) {
	// register struct stack_frame *fp asm("ebp");
	register struct stack_frame *fp asm("rbp");
	struct stack_frame *frame = fp;
	int i = 0;
	while(frame && i < maxframes) {
		abuf[i++] = frame->ret;
		frame = frame->next;
	}
	return i;
}
#endif

//
// simple stack trace implementaion - not for call from signal handler
//
static void PrintStackTraceFp(FILE *ofp, int start) {
	static void *addrlist[MAXFRAMES+1];
	int addrlen;
	char **symbollist;
	int i;

	myfputs("Stack Trace:\n", ofp);
	addrlen = backtrace(addrlist, MAXFRAMES);
	if(addrlen <= 3) {
		myfputs("Empty - may be corrupt\n", ofp);
		return;
	}
	/* I wish I could pass in the buffer */
	addrlen--;
	symbollist = backtrace_symbols(addrlist, addrlen);
	for(i = start; i < addrlen; i++) {
		myfputs(symbollist[i], ofp);
		fputc('\n', ofp);
	}
	return;
}

// fortran callable routine - note trailing _
void printstacktrace_(void) {
	PrintStackTraceFp(stdout, 2);
}

//
// generic signal handler for aborting signals
//
static void handle_sig(int signo, void *sip, my_sig_ucontext_t *scp) {
	FILE *ofp = stderr;
	int addrlen;
	char **symbollist;
	int i;

	myfputs("Caught signal ", stdout);
	fputc('0' + signo/10, stdout);
	fputc('0' + signo%10, stdout);
	fputc(' ', stdout);
	switch(signo) {
	case SIGSEGV:  
	    myfputs("SEGV", stdout); 
		/* printf("faulty address is %p, from %p ", ctx.cr2, ctx.eip); */
		break;
	case SIGBUS:   myfputs("BUS", stdout);  break;
	case SIGINT:   myfputs("INT", stdout);  break;
	case SIGABRT:  myfputs("ABRT", stdout);  break;
	case SIGUSR1:  myfputs("SIGUSR1", stdout); break;
	}
	myfputs("\nAttempting to print stack trace:\n", stdout);

	/* I'm adding it directly inside the signal handler so I use
	 * the undocumented arg ctx
	 */

	myfputs("Stack Trace:\n", ofp);
	addrlen = backtrace(addrlist, MAXFRAMES);
	if(addrlen <= 4) {
		myfputs("Empty - may be corrupt\n", ofp);
		myfputs("or May be missing -fno-omit-frame-pointer in build\n", ofp);
	} else {
#ifdef __APPLE__
		{
			//
			// For some reason, the stack in a signal handler lacks the usual
			// back pointer to the calling location.   See my_backtrace
			// above to see how it normally works.   The offical "backtrace"
			// function is broken and loses the innermost user execution point
			// These two lines of code restore it.
			//
			my_ucontext_t *uc = scp->uc_link;

			addrlist[2] = (void *)uc->regs.__rip;
		}
#endif
		/* I wish I could pass in the buffer */
		addrlen-=2; // trim off above MAIN 
		addrlen-=2; // trim off start two levels,  the signal handler and a special sigaction frame
		symbollist = backtrace_symbols(&addrlist[2], addrlen);
		for(i = 0; i < addrlen; i++) {
			myfputs(symbollist[i], ofp);
			fputc('\n', ofp);
		}
	}
	if(signo == SIGUSR1) {
		myfputs("Attempting to continue\n", ofp);
		return;
	}
	myfputs("Exiting...\n", stdout);
	exit(1);
}

void installsignalhandlers_(void) {
	int siglist[] = {SIGSEGV, SIGBUS, SIGABRT, SIGINT, SIGILL, SIGUSR1};
	int i;
	struct sigaction sa;

#if 0
   printf("Installing signal handlers SIGSEGV, SIGBUS, SIGABRT, SIGINT, SIGILL, SIGUSR1\n");
   printf("SIGUSR1 (30) continues after stack trace is printed\n");
#endif
	sa.sa_flags = 0;
	sigemptyset(&sa.sa_mask);
	sa.sa_handler = (void (*)) handle_sig;
	for(i = 0; i < sizeof(siglist)/sizeof(int); i++)  {
		sigaction(siglist[i], &sa, 0);
	}
}

/* Debugging routine to stop runaway processes */
void limitmemory_(int numgig) {
	struct rlimit rlim = {0};
	rlim_t max;
	max = numgig;
	max *= 1^30;
	rlim.rlim_max = rlim.rlim_cur = max;
	setrlimit(RLIMIT_DATA, &rlim);
}
