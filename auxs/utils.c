
#include "utils.h"
#include <float.h>
#include <unistd.h>
#include <signal.h>
#include <execinfo.h>

#ifdef __MACH__
#include <sys/time.h>
//clock_gettime is not implemented on OSX
int clock_gettime(int clk_id, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}
#endif

float maxs(size_t data[], size_t start, size_t n) {
	double max = DBL_MIN;
	for(int i = start; i<n; i++) {
		if (max < data[i])
			max = data[i];
	}
	return max;
}

float stdev(size_t data[], size_t n) {
	return stdevs(data, 0, n);
}

float stdevs(size_t data[], size_t start, size_t n) {
	float mean=0.0, sum_deviation=0.0;
	int i;
	for(i=start; i<n; ++i)
		mean+=data[i];
	mean=mean/(n-start);
	for(i=start; i<n; ++i)
		sum_deviation += (data[i]-mean)*(data[i]-mean);
	return sqrt(sum_deviation/(n-start));
}

double getvecvalue(const void *vec, const int n) {
	return ((double*)vec)[n];
}

double stdevd(double data[], size_t start, size_t n) {
	return stdevd_ex(&data, start, n, getvecvalue);
}

double stdevd_ex(void *data, size_t start, size_t n, double (*getv)(const void *, const int n)) {
	double mean=0.0, sum_deviation=0.0;
	int i;
	for(i=start; i<n; ++i)
		mean+=getv(data, i);
	mean=mean/(n-start);
	for(i=start; i<n; ++i)
		sum_deviation += (getv(data, i)-mean)*(getv(data, i)-mean);
	return sqrt(sum_deviation/(n-start));

}

void print_geojson_header() {
	fprintf(stderr, "{'type': 'FeatureCollection', 'features': [\n");
}

void print_geojson_mbr(const Envelope e, char *id) {
	fprintf(stderr, "{'type': 'Feature', 'geometry': {'type': 'Polygon', 'coordinates': [[");
	fprintf(stderr, "[%f, %f],", e.MinX, e.MinY);
	fprintf(stderr, "[%f, %f],", e.MaxX, e.MinY);
	fprintf(stderr, "[%f, %f],", e.MaxX, e.MaxY);
	fprintf(stderr, "[%f, %f],", e.MinX, e.MaxY);
	fprintf(stderr, "[%f, %f]",  e.MinX, e.MinY);
	fprintf(stderr, "]]}, 'properties': {'name': '%s'}},\n", id);
}

void print_geojson_footer() {
	fprintf(stderr, "]}\n");
}

int get_thread_num() {
	char *thread_num_str = getenv("THREAD_NUM");
	if (!thread_num_str)
		return 1;
	else
		return atoi(thread_num_str);
}

double runtime_diff_ms(struct timespec *start, struct timespec *end) {
	return ( end->tv_sec - start->tv_sec ) * 1000.0 + (double)( end->tv_nsec - start->tv_nsec ) / 1E6;
}

inline __attribute__((always_inline))
void print_progress_gauge(unsigned read, unsigned total) {
	const static char *progress_gauge_equal = "====================";
	const static char *progress_gauge_empty = "                    ";
	const static int size = 20;
	const static int resolution = 100/20;
	if ((read-1) % (total/100) != 0) return;
	int percent = ((read*100) / total);
	int pres = percent / resolution;
	fprintf(stderr, "  [%.*s%.*s] %3d%%\r", pres, progress_gauge_equal,
		size - pres, progress_gauge_empty, percent);
}

void handle_sigsegv(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stdout, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDOUT_FILENO);
  exit(1);
}

void init_geos() {
	initGEOS(geos_messages, geos_messages);
	signal(SIGSEGV, handle_sigsegv);
	signal(SIGABRT, handle_sigsegv);
}

GEOSContextHandle_t init_geos_context() {
	return initGEOS_r(geos_messages, geos_messages);
}
