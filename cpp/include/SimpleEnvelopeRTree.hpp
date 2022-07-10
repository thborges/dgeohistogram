/*
 * A cpp wrap for the rtree_lazy tree
 *
 *  Created on: 2022-07-10
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include "rtree-lazy.h"

template<class T, short m>
class SimpleEnvelopeRTree: public std::list<T> {
private:
    rtree_root *rtree;
public:
    SimpleEnvelopeRTree() {
		rtree = rtree_new_r0(m, 1);
    }

    ~SimpleEnvelopeRTree() {
        rtree_destroy(rtree);
    }

    void buildRTree() {
        for(const T& item: (*this))
            rtree_append(rtree, (GEOSGeometryH)&item, item.mbr);
    }

    std::list<T*> getIntersections(Envelope query) {
        std::list<T*> items;
		GList *results = rtree_window_search_mbr(rtree, query);
        GList *item;
        g_list_foreach(item, results) {
            rtree_leaf *leaf = (rtree_leaf *)item->data;
            items.push_back((T*)leaf->geo);
        }
        g_list_free(results);
        return items;
    }
};
