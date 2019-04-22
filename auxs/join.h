/*
 * join.h
 *
 *  Created on: 28/07/2014
 *      Author: thborges
 */

#ifndef JOIN_H_
#define JOIN_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "rtree.h"
#include "dataset.h"

/*Join two rtrees r and s, using RJ*/
dataset*		rtree_join_rj	(const rtree_root *r, const rtree_root *s, rtree_join_stat *stats);

/*Join a dataset r and an rtree s, using INL*/
dataset*		rtree_join_inl	(dataset *r, const rtree_root *s, rtree_join_stat *stats, enum JoinPredicateCheck pcheck);

/*Join two datasets r and s, using HJ*/
dataset*		rtree_join_hj	(dataset *r, dataset *s, rtree_join_stat *stats, enum JoinPredicateCheck pcheck);

/*Parallel Join two datasets r and s, using HJ*/
extern dataset*	rtree_pjoin_hj	(dataset *r, dataset *s, rtree_join_stat *stats, enum JoinPredicateCheck pcheck);

/*Parallel Join two rtrees r and s, using RJ*/
extern dataset*	rtree_pjoin_rj	(const rtree_root *r, const rtree_root *s, rtree_join_stat *stats);

/*Parallel Join a dataset r and an rtree s, using INL*/
extern dataset* rtree_pjoin_inl	(dataset* r, const rtree_root *s, rtree_join_stat *stats, enum JoinPredicateCheck pcheck);

#ifdef __cplusplus
}
#endif

#endif /* JOIN_H_ */
