/*
 * A deque with binary search for point object
 *
 *  Created on: 2022-07-09
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include "rtree-lazy.h"
#include <deque>
#include <algorithm>

template<
    class T,
    class Allocator = std::allocator<T>>
class DequeBinarySearchPoint: public std::deque<T, Allocator> {
public:
    void sortForSL() {
        std::sort(this->begin(), this->end(), [](const T & a, const T & b) -> bool {
            return a.x < b.x;
        });
    }

    std::list<T*> getIntersectionsSL(const Envelope &query) {
        std::list<T*> items;
        if (this->size() == 0)
            return items;

        // binary search for the first point after query.MinX
        long first = 0;
        long last = this->size()-1;
        long middle;
        while (first <= last) {
            middle = first+(last-first)/2;
            const T& current = (*this)[middle];
            if (current.x == query.MinX) {
                first = middle;
                break;
            } else if (current.x < query.MinX)
                first = middle+1;
            else if (current.x > query.MinX)
                last = middle-1;
        }
        // binary search doesn't ensure that the element at 'first' is
        // the first with x equal to query.MinX
        while (first > 0 && (*this)[first].x >= query.MinX)
            first--;

        for(int i = first; i < this->size(); i++) {
            T* current = &(*this)[i];
            if (current->x > query.MaxX)
                break;
            if (query.contains(current->x, current->y))
                items.push_back(current);
        }
        return items;
    }
};
