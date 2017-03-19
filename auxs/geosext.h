/*
 * geosext.h
 *
 *  Created on: Jun 18, 2014
 *      Author: thborges
 */

#ifndef GEOSEXT_H_
#define GEOSEXT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <geos_c.h>

#ifndef GEOSIEnvelope
typedef struct GEOSEnvelope_t GEOSIEnvelope;
#endif

extern const GEOSIEnvelope *GEOSInternalEnvelope(const GEOSGeometry *g1);
extern void GEOSEnvelopeGetXY(const GEOSGeometry *g1, double *x1, double *x2, double *y1, double *y2);
extern size_t GEOSGetNumPoints(const GEOSGeometry *e);
extern void geos_messages(const char *fmt, ...);

extern void *GEOSCreateWKBBuffer();
extern void *GEOSCreateWKTBuffer();
extern const char *GEOSGeomToWKB(const GEOSWKBWriter *w, const void *buff, const GEOSGeometry *g, int *size);
extern const char *DgeoGeomToWKT(const GEOSWKTWriter *w, const void *buff, const GEOSGeometry *g);
extern GEOSWKBWriter *GEOSCreateWKBWritter(GEOSContextHandle_t extHandle);
extern GEOSWKTWriter *GEOSCreateWKTWritter();
extern void GEOSDestroyWKBWritter(GEOSWKBWriter *w, void *buff);
extern void GEOSDestroyWKTWritter(GEOSWKTWriter *w, void *buff);

#ifdef __cplusplus
}
#endif

#endif /* GEOSEXT_H_ */
