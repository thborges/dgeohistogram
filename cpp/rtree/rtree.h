#ifndef RTREE_H
#define RTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <geos_c.h>
#include <stdbool.h>
#include <assert.h>
#include <sys/queue.h>
#include "utils.h"
#include "glibwrap.h"
#include "ogrext.h"

#define VPOINTER(x) &x[0]

#define minthreshold 0.4
#define reinsertthreshold 0.3
#define minm(M) (int)ceil(M*minthreshold)
#define maxm(M) (int)ceil(M*(1-minthreshold))+1

#define SWAP(x, y, T) do { T temp##x##y = x; x = y; y = temp##x##y; } while (0)

enum NodeTypeEnum {DIRECTORY, LEAF};

struct leaf {
    bool cloned;
    long long gid;
    unsigned points;
    EnvelopeC mbr;
    //GEOSGeometryH chull;
    GEOSGeometryH geo;
    const GEOSPreparedGeometry *pgeo;
};
typedef struct leaf rtree_leaf;

struct rtree_node {
	EnvelopeC mbr;

	union {
		struct rtree_node **dirs;
		rtree_leaf *leaves;
	};

	int index;
	int server;

	enum NodeTypeEnum type;
	short used;
};

typedef struct rtree_node rtree_node;


typedef struct rtree_root_ {
	int m;
	int servers;
	int rr_next_server;
	int levels;
	rtree_node *root;
	int (*choose_subtree_func)(struct rtree_root_ *, rtree_node *, const EnvelopeC e);
	rtree_node *(*rtree_split_func)(struct rtree_root_ *, rtree_node *, rtree_node *, rtree_node *, int);
	rtree_node *(*rtree_split_leaf_func)(struct rtree_root_ *, rtree_node *, rtree_node *, const GEOSGeometryH, const EnvelopeC);
	void (*reinsert_on_parent)(rtree_node *, struct rtree_root_ *, GList **, const GEOSGeometryH, const EnvelopeC);
	void (*rtree_redivide)(rtree_node *, rtree_node *, const int);
	void (*rtree_redivide_leaf)(rtree_node *, rtree_node *, const int);
} rtree_root;

typedef struct {
	GEOSGeometryH geo;
	const GEOSPreparedGeometry *pgeo;
	GEOSGeometryH chull;
	EnvelopeC mbr;
} rtree_window;

typedef struct {
	float trueleaves;
	float falseleaves;
	int truedirs;
	int falsedirs;
	int rmessages;
	int lmessages;
	int geomchecked;
} rtree_window_stat;

typedef struct {
	size_t io_pages;
	double cache_miss;
	unsigned mbrs_compared;
	unsigned geos_compared;
	unsigned mbrs_compared_false;
	unsigned geos_compared_false;
	int queue1_size;
	int queue2_size;
	int lines;
} rtree_join_stat;

typedef struct {
	rtree_node *nr;
	rtree_node *ns;
} rtree_node_pair;

typedef struct {
    int leaves;
    float subleaves;
    int dirs;
    int items;
} count_node_struct;

/*Create a new empty rtree index with fanout M.*/
rtree_root *rtree_new(
	const int m,
	const int servers, 
	int         (*choose_subtree_func)   (rtree_root *, rtree_node *, const EnvelopeC e),
	rtree_node *(*rtree_split_func)      (rtree_root *, rtree_node *, rtree_node *, rtree_node *, int),
	rtree_node *(*rtree_split_leaf_func) (struct rtree_root_ *, rtree_node *, rtree_node *, const GEOSGeometryH, const EnvelopeC e),
	void        (*reinsert_on_parent)    (rtree_node *, struct rtree_root_ *, GList **, const GEOSGeometryH, const EnvelopeC),
	void        (*rtree_redivide)        (rtree_node *, rtree_node *, const int),
	void        (*rtree_redivide_leaf)   (rtree_node *, rtree_node *, const int)

);

void		rtree_destroy			(rtree_root *r);

/*Add a new geometry to the rtree index r.*/
void		rtree_append			(rtree_root *r, const GEOSGeometryH g, const EnvelopeC mbr);

/*Return the height of the rtree r.*/
int 		rtree_height			(const rtree_root *r);

/*Find geometries on rtree r, wich intercepts Envelope window.*/
GList		*rtree_window_search	(const rtree_root *r, const GEOSGeometryH window, const EnvelopeC env, rtree_window_stat *stats);

/*Find candidate nodes on rtree r, wich intercepts Envelope window.*/
GList		*rtree_window_csearch	(const rtree_root *r, const GEOSGeometryH window, const EnvelopeC env, rtree_window_stat *stats);

/*Compute the total mbr of the node n*/
EnvelopeC	rtree_compute_mbr		(const rtree_node *n) G_GNUC_WARN_UNUSED_RESULT;

/* auxiliary internal functions */
rtree_node *rtree_new_node(rtree_root *root, const enum NodeTypeEnum type);
void rtree_append_internal(rtree_root *r, const GEOSGeometryH g, const EnvelopeC gmbr, rtree_node *at_reinsert);

/* page buffer */
//extern lrubuffer *lru;

#ifdef __cplusplus
}
#endif

#endif
