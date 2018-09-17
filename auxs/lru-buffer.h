
#ifndef LRU_BUFFER
#define LRU_BUFFER

#include "glibwrap.h"
#include "uthash.h"

typedef struct {
	gpointer id;
	/*char *data;*/
	UT_hash_handle hh;
} lru_page;

typedef struct {
	lru_page *hasht;
	GList *list;
	GList *list_last;
	size_t list_size;
	size_t pages;
	size_t get_extmem;
	size_t get_buffer;
} lrubuffer;

lrubuffer *lrubuffer_new(size_t pages);
void lrubuffer_destroy(lrubuffer *lb);
void *lrubuffer_get(lrubuffer *lb, gpointer page);
#endif

