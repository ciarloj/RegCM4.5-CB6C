#include "private.h"

#if ( defined THREADED_OMP )

#include <omp.h>

/*
** threadinit: Initialize locking capability and set number of threads
**
** Output arguments:
**   nthreads:   number of threads
**   maxthreads: number of threads (these don't differ under OpenMP)
*/

int threadinit (int *nthreads, int *maxthreads)
{
  /* In OMP case, maxthreads and nthreads are the same number */

  *maxthreads = omp_get_max_threads ();
  *nthreads = *maxthreads;
  if (get_thread_num (nthreads, maxthreads) > 0)
    return GPTLerror ("GPTLthreadinit: MUST be called only by master thread");
  return 0;
}

/*
** threadfinalize: no-op under OpenMP
*/

void threadfinalize ()
{
}

/*
** get_thread_num: determine thread number of the calling thread
**
** Input args:
**   nthreads:   number of threads
**   maxthreads: number of threads (unused in OpenMP case)
**
** Return value: thread number (success) or GPTLerror (failure)
*/

int get_thread_num (int *nthreads, int *maxthreads)
{
  int t;       /* thread number */

  if ((t = omp_get_thread_num ()) >= *nthreads)
    return GPTLerror ("get_thread_num: returned id %d exceed numthreads %d\n",
		      t, *nthreads);

  return t;
}

#elif ( defined THREADED_PTHREADS )

#define MAX_THREADS 128

#include <pthread.h>
#include <stdlib.h>
static int lock_mutex (void);      /* lock a mutex for entry into a critical region */
static int unlock_mutex (void);    /* unlock a mutex for exit from a critical region */

static pthread_mutex_t t_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_t *threadid;

/*
** threadinit: Set number of threads and max number of threads
**
** Output arguments:
**   nthreads:   number of threads
**   maxthreads: number of threads (these don't differ under OpenMP)
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int threadinit (int *nthreads, int *maxthreads)
{
  int nbytes;

  /* Manage the threadid array which maps physical thread id's to logical id's */

  nbytes = MAX_THREADS * sizeof (pthread_t);
  if ( ! (threadid = (pthread_t *) malloc (nbytes)))
    return GPTLerror ("threadinit: malloc failure for %d items\n", MAX_THREADS);

  /*
  ** Initialize nthreads to 1 and define the threadid array now that initialization 
  ** is done. The actual value will be determined as get_thread_num is called.
  */

  threadid[0] = pthread_self ();
  *nthreads = 1;
  *maxthreads = MAX_THREADS;

  return 0;
}

/*
** threadfinalize: free allocated space
*/

void threadfinalize ()
{
  free (threadid);
}

/*
** get_thread_num: determine zero-based thread number of the calling thread
**
** Input args:
**   nthreads:   number of threads
**
** Input/output args:
**   maxthreads: max number of threads
**
** Return value: thread number (success) or GPTLerror (failure)
*/

int get_thread_num (int *nthreads, int *maxthreads)
{
  int n;                 /* return value: loop index over number of threads */
  pthread_t mythreadid;  /* thread id from pthreads library */

  mythreadid = pthread_self ();

  if (lock_mutex () < 0)
    return GPTLerror ("get_thread_num: mutex lock failure\n");

  /*
  ** Loop over known physical thread id's.  When my id is found, map it 
  ** to logical thread id for indexing.  If not found return a negative 
  ** number.
  ** A critical region is necessary because acess to
  ** the array threadid must be by only one thread at a time.
  */

  for (n = 0; n < *nthreads; n++)
    if (pthread_equal (mythreadid, threadid[n]))
      break;

  /*
  ** If our thread id is not in the known list, add to it after checking that
  ** we do not have too many threads.
  */

  if (n == *nthreads) {
    if (*nthreads >= MAX_THREADS) {
      if (unlock_mutex () < 0)
	printf ("get_thread_num: mutex unlock failure\n");
      return GPTLerror ("get_thread_num: nthreads=%d is too big Recompile "
			"with larger value of MAX_THREADS\n", *nthreads);
    }    
    threadid[n] = mythreadid;
    ++*nthreads;
  }
    
  if (unlock_mutex () < 0)
    return GPTLerror ("get_thread_num: mutex unlock failure\n");

  return n;
}

/*
** lock_mutex: lock a mutex for private access
*/

static int lock_mutex ()
{
  if (pthread_mutex_lock (&t_mutex) != 0)
    return GPTLerror ("pthread_lock_mutex failure\n");
  return 0;
}

/*
** unlock_mutex: unlock a mutex from private access
*/

static int unlock_mutex ()
{
  if (pthread_mutex_unlock (&t_mutex) != 0)
    return GPTLerror ("pthread_unlock_mutex failure\n");
  return 0;
}

#else

/*
** Unthreaded case
*/

int threadinit (int *nthreads, int *maxthreads)
{
  *nthreads = 1;
  *maxthreads = 1;
  return 0;
}

int get_thread_num (int *nthreads, int *maxthreads)
{
  return 0;
}
void threadfinalize ()
{
}
#endif
