/*
 * glibwrap.c
 *
 *  Created on: 14/07/2015
 *      Author: thborges
 */

#include "glibwrap.h"

void g_free(void *pointer) {
	free(pointer);
}

void *g_memdup(const void *mem, unsigned int byte_size) {
	void *result = g_new(unsigned char, byte_size);
	return memcpy(result, mem, byte_size);
}

GList *g_new_list(void *data) {
	GList *lst = g_new0(GList, 1);
	lst->data = data;
	return lst;
}

GList *g_list_append(GList *lst, void *data) {
	if (lst == NULL) 
		return g_new_list(data);

	GList *last = lst;
	while(last->next)
		last = last->next;

	GList *newi = g_new(GList, 1);
	newi->data = data;
	newi->next = NULL;
	newi->priour = last;

	last->next = newi;

	return lst;
}

GList *g_list_concat(GList *lst, GList *lst2) {
	GList *last = lst;
	while(last->next)
		last = last->next;
	last->next = lst2;
	lst2->priour = last;
	return last;
}

GList *g_list_prepend(GList *lst, void *data) {
	GList *newi = g_new(GList, 1);
	newi->next = lst;
	newi->priour = NULL;
	newi->data = data;
	if (lst)
		lst->priour = newi;
	return newi;
}

GList *g_list_find(GList *lst, void *data) {
	while(lst) {
		if (lst->data == data) return lst;
		lst = lst->next;
	}
	return NULL;
}

GList *g_list_remove(GList *lst, void *data) {
	GList *link = g_list_find(lst, data);
	if (link) {
		if (link->next)
			link->next->priour = link->priour;
		if (link->priour)
			link->priour->next = link->next;
		else
			lst = link->next;
	}
	return lst;
}

GList *g_list_remove_link(GList *lst, GList *link) {
	if (link->next)
		link->next->priour = link->priour;
	if (link->priour)
		link->priour->next = link->next;
	else
		lst = link->next;
	return lst;
}

GList *g_list_delete_link(GList *lst, GList *link) {
	lst = g_list_remove_link(lst, link);
	g_free(link);
	return lst;
}

void g_list_free_full(GList *lst, void (*destroy)(void *data)) {
	while (lst) {
		GList *aux = lst;
		lst = lst->next;
		destroy(aux->data);
		g_free(aux);
	}
}

void g_list_free(GList *lst) {
	while (lst) {
		GList *aux = lst;
		lst = lst->next;
		g_free(aux);
	}
}

unsigned g_list_length(GList *lst) {
	unsigned result = 0;
	while (lst) {
		result++;
		lst = lst->next;
	}
	return result;
}
