
#ifndef UTILS_H
#define UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ogrext.h>

//#define INT_RAND_MINMAX(min, max) (((double) random() / (RAND_MAX+1)) * (max-min+1) + min)
#define INT_RAND(max) (rand() % max + 1)

#define	STAILQ_HEADS(name, type)                                \
struct name {                                                   \
    struct type *stqh_first;/* first element */                 \
    struct type **stqh_last;/* addr of last next element */		\
    size_t count;/* number of elements */                       \
}

#define	STAILQ_INITS(head) STAILQ_INIT(head); (head)->count = 0

#define STAILQ_INSERT_TAILS(head, elm, field) STAILQ_INSERT_TAIL(head, elm, field); (head)->count++ 

#define STAILQ_CONCATS(head1, head2) do { \
	(head1)->count += (head2)->count;     \
	STAILQ_CONCAT(head1, head2);          \
  } while (0)

#define STAILQ_FREE_FULL(type, head, field) do {         \
    type *lp__ = STAILQ_FIRST(head);                    \
    while (lp__ != NULL) {                              \
        type *lpa__ = STAILQ_NEXT(lp__, field);          \
        g_free(lp__);                                   \
        lp__ = lpa__;                                   \
    }} while (0)

#define STAILQ_REMOVE_HEADS(head, field) (head)->count--; STAILQ_REMOVE_HEAD(head, field)

float stdev(size_t data[], size_t n);
float stdevs(size_t data[], size_t start, size_t n);
double stdevd(double data[], size_t start, size_t n);
double stdevd_ex(void *data, size_t start, size_t n, double (*getv)(const void *, const int n));
float maxs(size_t data[], size_t start, size_t n);
void print_geojson_mbr(const Envelope e, char *id);
void print_geojson_footer();
void print_geojson_header(); 
int get_thread_num();

#ifdef __MACH__
#include <sys/time.h>
//clock_gettime is not implemented on OSX
int clock_gettime(int /*clk_id*/, struct timespec* t);
#endif
#define CLOCK_REALTIME 0
#define CLOCK_PROCESS_CPUTIME_ID 2

double runtime_diff_ms(struct timespec *start, struct timespec *end);
void print_progress_gauge(unsigned read, unsigned total);
void geos_messages(const char *fmt, ...);
void handle_sigsegv(int sig);
void init_geos();
GEOSContextHandle_t init_geos_context();

// declare and call qsort_r on OSX and Linux, in a portable way
#ifdef __MACH__
#define decl_qsort_p_cmp(fname, x, y, thunk) int fname(void *thunk, const void *x, const void *y)
#define qsort_p(base, nmemb, size, compar, thunk) qsort_r(base, nmemb, size, thunk, compar)
#elif __linux__
#define decl_qsort_p_cmp(fname, x, y, thunk) int fname(const void *x, const void *y, void *thunk)
#define qsort_p(base, nmemb, size, compar, thunk) qsort_r(base, nmemb, size, compar, thunk)
#endif

#ifdef __cplusplus
}
#endif

#endif

