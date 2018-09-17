#include "lru-buffer.h"

lrubuffer *lrubuffer_new(size_t pages) {
	lrubuffer *lb = g_new(lrubuffer, 1);
	lb->hasht = NULL;
	lb->list = NULL;
	lb->list_last = NULL;
	lb->list_size = 0;
	lb->pages = pages;
	lb->get_extmem = 0;
	lb->get_buffer = 0;
	return lb;
}

void lrubuffer_destroy(lrubuffer *lb) {
	lru_page *current_page, *tmp;
	HASH_ITER(hh, lb->hasht, current_page, tmp) {
		HASH_DEL(lb->hasht, current_page);
	}
	g_list_free_full(lb->list, g_free);
	g_free(lb->hasht);
	g_free(lb);
	lb = NULL;
}

void *lrubuffer_get(lrubuffer *lb, gpointer page) {
	return NULL;

	lru_page *p;
	HASH_FIND_PTR(lb->hasht, &page, p);
	if (p != NULL) {
		GList *page = g_list_find(lb->list, p); //not so efficient. for test only
		if (page == lb->list_last)
			lb->list_last = g_list_previous(lb->list_last);

		lb->list = g_list_remove_link(lb->list, page);
		lb->list = g_list_concat(page, lb->list);
		lb->get_buffer++; // get from internal buffer
	}
	else {
		lb->get_extmem++; // get from external memory
		lru_page *new_page = NULL;

		// full buffer. Remove one page.
		if (lb->list_size >= lb->pages) {
			GList *tail = lb->list_last;
			lru_page *tailp = tail->data;

			lb->list_last = g_list_previous(lb->list_last);
			lb->list = g_list_delete_link(lb->list, tail);
			
			HASH_FIND_PTR(lb->hasht, &tailp->id, p);
			if (p != NULL) {
				HASH_DEL(lb->hasht, p);
			}
			new_page = tailp;
		}
		else {
			new_page = g_new(lru_page, 1);
			lb->list_size++;
		}
		new_page->id = page;

		HASH_ADD_PTR(lb->hasht, id, new_page);
		lb->list = g_list_prepend(lb->list, new_page);
		if (lb->list_last == NULL)
			lb->list_last = lb->list;
	}
	return page; //fake
}

