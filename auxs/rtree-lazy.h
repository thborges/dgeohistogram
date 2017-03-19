
#ifndef RTREE_LAZY_H
#define RTREE_LAZY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "rtree.h"
#include "ogrext.h"
#include <assert.h>

// R0-tree definitions
#define ENV_WIDTH(r) (r.MaxX - r.MinX)
#define ENV_HEIGHT(r) (r.MaxY - r.MinY)
#define R0_QUALITY(r) ((1.0/(ENV_WIDTH(r)*ENV_HEIGHT(r))) * powf( (MIN(ENV_WIDTH(r),ENV_HEIGHT(r))/MAX(ENV_WIDTH(r),ENV_HEIGHT(r))), 0.5))
#define R0_GAIN(r1, r2) (1-(R0_QUALITY(r1)/R0_QUALITY(r2)))
#define R0_LOSS(r1, r2) R0_GAIN(r2,r1)

#define PROP_MARGIN(r) (MIN(ENV_WIDTH(r),ENV_HEIGHT(r))/MAX(ENV_WIDTH(r),ENV_HEIGHT(r)))
#define MIN_PERIMETER(r1, r2) (MIN(ENVELOPE_PERIMETER(r1), ENVELOPE_PERIMETER(r2)))
#define MAX_PERIMETER(r1, r2) (MAX(ENVELOPE_PERIMETER(r1), ENVELOPE_PERIMETER(r2)))

//#define SPLIT_QUALITY(overlap, area, r1, r2) ((1.0/((overlap == 0.0 ? 1.0e-10 : overlap) / area)))
#define SPLIT_QUALITY(overlap, area, r1, r2) ((1.0/((overlap == 0.0 ? 1.0e-10 : overlap) * area)))

// * powf( MIN_PERIMETER(r1,r2) / MAX_PERIMETER(r1,r2), 0.5)))
//#define SPLIT_QUALITY(overlap, area, r1, r2) (1.0/(area) * (overlap == 0.0 ? 1 : powf(1.0/(overlap), 0.9)))// * powf( MIN_PERIMETER(r1,r2) / MAX_PERIMETER(r1,r2), 0.05))

#define R1_QUALITY(r, over) ((1.0/(ENV_WIDTH(r)*ENV_HEIGHT(r))) * (overlap == 0.0 ? 1 : powf( 1.0/(over), 0.9)))
// * powf( (MIN(ENV_WIDTH(r),ENV_HEIGHT(r))/MAX(ENV_WIDTH(r),ENV_HEIGHT(r))), 0.1) )
#define R1_GAIN(r1, overr1, r2, overr2) (1-(R1_QUALITY(r1, overr1)/R1_QUALITY(r2, overr2)))
#define R1_LOSS(r1, overr1, r2, overr2) R1_GAIN(r2, overr2, r1, overr1)

rtree_root *rtree_new_r0(const int m, const int servers);
rtree_root *rtree_new_rlazy(const int m, const int servers);

#ifdef __cplusplus
}
#endif

#endif
