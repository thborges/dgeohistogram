#ifndef OGREXT_H
#define OGREXT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <geos_c.h>
#include <ogr_api.h>
#include <emmintrin.h>
//#include <iacaMarks.h>
#include <assert.h>

struct Envelope { 
	union {
		struct {
			__m128d Min;
			__m128d Max;
		};
		struct {
			double MinX;
			double MinY;
			double MaxX;
			double MaxY;
		};
	};
};

typedef struct Envelope Envelope;

extern Envelope emptymbr;

/*Returns the intersection envelope of r an s*/
const Envelope EnvelopeIntersection(const Envelope r, const Envelope s);
const Envelope EnvelopeIntersection2(const Envelope r, const Envelope s);

/*Convert a Envelope to a OGRGeometry*/
OGRGeometryH EnvelopeToGeometry(const Envelope r);

void GetEnvelopeCoordinates(Envelope p, double *MinX, double *MinY, double *MaxX, double *MaxY);
Envelope GetEnvelope(double MinX, double MinY, double MaxX, double MaxY);

void GEOSGetEnvelope(const GEOSGeometry *geo, Envelope *env);

GEOSGeometry *OGR_G_ExportToGEOS( const OGRGeometryH ogrg);

Envelope OGRGetEnvelope(OGRGeometryH ogrg);
int OGRGetNumPoints(OGRGeometryH ogrg);

#define ENVELOPE_MERGE(mbr, e) do { \
mbr.MinX = MIN(mbr.MinX, e.MinX); \
mbr.MinY = MIN(mbr.MinY, e.MinY); \
mbr.MaxX = MAX(mbr.MaxX, e.MaxX); \
mbr.MaxY = MAX(mbr.MaxY, e.MaxY); } while(0)

#define ENVELOPE_BUFFER(mbr, sizew, sizeh) \
mbr.MinX -= sizew; \
mbr.MinY -= sizeh; \
mbr.MaxX += sizew; \
mbr.MaxY += sizeh

#define ENVELOPE_INIT(mbr) mbr = emptymbr

#define ENVELOPE_AREA(mbr) (((mbr).MaxX - (mbr).MinX) * ((mbr).MaxY - (mbr).MinY))

#define ENVELOPE_PERIMETER(mbr) (2*(mbr.MaxX - mbr.MinX) + 2*(mbr.MaxY - mbr.MinY))

//#define ENVELOPE_INTERSECTS(mbr1, mbr2) (!(mbr2.MinX > mbr1.MaxX || mbr2.MaxX < mbr1.MinX || mbr2.MinY > mbr1.MaxY || mbr2.MaxY < mbr1.MinY))
//#define ENVELOPE_INTERSECTS(mbr1, mbr2) (mbr1.MinX <= mbr2.MaxX && mbr1.MaxX >= mbr2.MinX && mbr1.MinY <= mbr2.MaxY && mbr1.MaxY >= mbr2.MinY)
//#define ENVELOPE_INTERSECTS(mbr1, mbr2) (mbr1.MinX <= mbr2.MaxX && mbr1.MinY <= mbr2.MaxY && mbr2.MinX <= mbr1.MaxX && mbr2.MinY <= mbr1.MaxY)

/*static inline __attribute__((always_inline))
char ENVELOPE_INTERSECTS(const Envelope mbr1, const Envelope mbr2) {
	//IACA_START
	char r = (!(mbr2.MinX > mbr1.MaxX || mbr2.MaxX < mbr1.MinX || mbr2.MinY > mbr1.MaxY || mbr2.MaxY < mbr1.MinY));
	//char r = (mbr1.MinX <= mbr2.MaxX && mbr1.MinY <= mbr2.MaxY && mbr2.MinX <= mbr1.MaxX && mbr2.MinY <= mbr1.MaxY);
	//IACA_END
	return r;
}*/

void print_geojson_mbr_local(const Envelope e, const char *id);

//#define ENVELOPE_INTERSECTS(mbr1, mbr2) (mbr1.MinX <= mbr2.MaxX && mbr1.MinY <= mbr2.MaxY && mbr2.MinX <= mbr1.MaxX && mbr2.MinY <= mbr1.MaxY)
static inline __attribute__((always_inline))
char ENVELOPE_INTERSECTS(const Envelope mbr1, const Envelope mbr2) {
	//IACA_START
	char result;
	/*__asm__( // using or
		"cmpnlepd %3, %4;" 
		"movmskpd %4, %%ecx;"
		"cmp $0x0, %%cl;"
		"jnz 2f;"
		"cmpnlepd %2, %1;"
		"movmskpd %1, %%ecx;"
		"cmp $0x0, %%cl;"
		"2: setz %0;"
		: "=r" (result)
		: "x" (mbr1.Min), "x" (mbr2.Max), "x" (mbr1.Max), "x" (mbr2.Min)
		: "%ecx");*/

	/*__asm__( // using and
		"cmplepd %3, %4;" 
		"movmskpd %4, %%ecx;"
		"cmp $0x3, %%cl;"
		"jne 2f;"
		"cmplepd %2, %1;"
		"movmskpd %1, %%ecx;"
		"cmp $0x3, %%cl;"
		"2: setz %0;"
		: "=r" (result)
		: "x" (mbr1.Min), "x" (mbr2.Max), "x" (mbr1.Max), "x" (mbr2.Min)
		: "%ecx");*/

	/*__asm__ ( // without branch
		"cmplepd %2, %1;"
		"cmplepd %3, %4;" 
		"andpd %1, %4;"
		"movmskpd %4, %%ebx;"
		"cmp $0x3, %%bl;"
		"setz %0;"
		: "=r" (result)
		: "x" (mbr1.Min), "x" (mbr2.Max), "x" (mbr1.Max), "x" (mbr2.Min)
		: "%ebx");*/
	
	result = (mbr1.MinX <= mbr2.MaxX && mbr1.MinY <= mbr2.MaxY && mbr2.MinX <= mbr1.MaxX && mbr2.MinY <= mbr1.MaxY);
	/*if (r != result) {
		print_geojson_mbr_local(mbr1, "0");
		print_geojson_mbr_local(mbr2, "1");
		printf("\n");
	}*/
	//IACA_END
	return result;
}

static inline __attribute__((always_inline))
char ENVELOPE_INTERSECTS_EULERSKEW(const Envelope mbr1, const Envelope mbr2) {
	//IACA_START
	char result;
	/*__asm__( // using or
		"cmpnlepd %3, %4;" 
		"movmskpd %4, %%ecx;"
		"cmp $0x0, %%cl;"
		"jnz 2f;"
		"cmpnlepd %2, %1;"
		"movmskpd %1, %%ecx;"
		"cmp $0x0, %%cl;"
		"2: setz %0;"
		: "=r" (result)
		: "x" (mbr1.Min), "x" (mbr2.Max), "x" (mbr1.Max), "x" (mbr2.Min)
		: "%ecx");*/

	/*__asm__( // using and
		"cmplepd %3, %4;" 
		"movmskpd %4, %%ecx;"
		"cmp $0x3, %%cl;"
		"jne 2f;"
		"cmplepd %2, %1;"
		"movmskpd %1, %%ecx;"
		"cmp $0x3, %%cl;"
		"2: setz %0;"
		: "=r" (result)
		: "x" (mbr1.Min), "x" (mbr2.Max), "x" (mbr1.Max), "x" (mbr2.Min)
		: "%ecx");*/

	/*__asm__ ( // without branch
		"cmplepd %2, %1;"
		"cmplepd %3, %4;" 
		"andpd %1, %4;"
		"movmskpd %4, %%ebx;"
		"cmp $0x3, %%bl;"
		"setz %0;"
		: "=r" (result)
		: "x" (mbr1.Min), "x" (mbr2.Max), "x" (mbr1.Max), "x" (mbr2.Min)
		: "%ebx");*/
	
	result = (mbr1.MinX < mbr2.MaxX && mbr1.MinY < mbr2.MaxY && mbr2.MinX < mbr1.MaxX && mbr2.MinY < mbr1.MaxY);
	/*if (r != result) {
		print_geojson_mbr_local(mbr1, "0");
		print_geojson_mbr_local(mbr2, "1");
		printf("\n");
	}*/
	//IACA_END
	return result;
}

#define ENVELOPE_CONTAINS(o, i) (i.MinX >= o.MinX && i.MaxX <= o.MaxX && i.MinY >= o.MinY && i.MaxY <= o.MaxY)
#define ENVELOPE_CONTAINS_EULERSKEW(o, i) (i.MinX > o.MinX && i.MaxX < o.MaxX && i.MinY > o.MinY && i.MaxY < o.MaxY)
#define ENVELOPE_CONTAINSP(e, x, y) (x >= e.MinX && x <= e.MaxX && y >= e.MinY && y <= e.MaxY)
#define ENVELOPE_CONTAINSP2(e, x, y) (x > e.MinX && x < e.MaxX && y > e.MinY && y < e.MaxY)


#ifdef __cplusplus
}
#endif

#endif

