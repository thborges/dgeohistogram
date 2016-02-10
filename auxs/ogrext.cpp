//#include <ogr_api.h>
#include <ogr_geometry.h>
#include <geos_c.h>
#include "ogrext.h"
#include "geosext.h"

Envelope emptymbr = {0.0, 0.0, 0.0, 0.0};

extern "C" const Envelope EnvelopeIntersection(const Envelope r, const Envelope s);
extern "C" const Envelope EnvelopeIntersection2(const Envelope r, const Envelope s);
extern "C" void GEOSGetEnvelope(const GEOSGeometry *geo, Envelope* env);
extern "C" GEOSGeometry *OGR_G_ExportToGEOS( const OGRGeometryH ogrg);
extern "C" Envelope OGRGetEnvelope(OGRGeometryH ogrg);
extern "C" int OGRGetNumPoints(OGRGeometryH ogrg);
extern "C" void GetEnvelopeCoordinates(Envelope p, double *MinX, double *MinY, double *MaxX, double *MaxY);
extern "C" Envelope GetEnvelope(double MinX, double MinY, double MaxX, double MaxY);

const Envelope EnvelopeIntersection(const Envelope r, const Envelope s) {
    if (ENVELOPE_INTERSECTS(r, s)) {
        Envelope result;
        result.MinX = r.MinX > s.MinX ? r.MinX : s.MinX;
        result.MinY = r.MinY > s.MinY ? r.MinY : s.MinY;
        result.MaxX = r.MaxX < s.MaxX ? r.MaxX : s.MaxX;
        result.MaxY = r.MaxY < s.MaxY ? r.MaxY : s.MaxY;
        return result;
    }
    else
        return emptymbr;
}

/*
 * Don't check intersection between r and s. Need to be done on calling code 
 */
const Envelope EnvelopeIntersection2(const Envelope r, const Envelope s) {
    Envelope result = {
    		r.MinX > s.MinX ? r.MinX : s.MinX,
    		r.MinY > s.MinY ? r.MinY : s.MinY,
    		r.MaxX < s.MaxX ? r.MaxX : s.MaxX,
    		r.MaxY < s.MaxY ? r.MaxY : s.MaxY};
    return result;
}

void GEOSGetEnvelope(const GEOSGeometry *geo, Envelope* env) {
	GEOSEnvelopeGetXY(geo, &env->MinX, &env->MaxX, &env->MinY, &env->MaxY);
}

GEOSGeometry *OGR_G_ExportToGEOS( const OGRGeometryH ogrg) {
	return ((OGRGeometry *) ogrg)->exportToGEOS(NULL);
}

void print_geojson_mbr_local(const Envelope e, const char *id) {
	fprintf(stderr, "{'type': 'Feature', 'geometry': {'type': 'Polygon', 'coordinates': [[");
	fprintf(stderr, "[%f, %f],", e.MinX, e.MinY);
	fprintf(stderr, "[%f, %f],", e.MaxX, e.MinY);
	fprintf(stderr, "[%f, %f],", e.MaxX, e.MaxY);
	fprintf(stderr, "[%f, %f],", e.MinX, e.MaxY);
	fprintf(stderr, "[%f, %f]",  e.MinX, e.MinY);
	fprintf(stderr, "]]}, 'properties': {'name': '%s'}},\n", id);
}

Envelope OGRGetEnvelope(OGRGeometryH ogrg) {
	Envelope result;
	OGREnvelope ev;
	OGR_G_GetEnvelope(ogrg, &ev);
	return Envelope{ev.MinX, ev.MinY, ev.MaxX, ev.MaxY};
}

int OGRGetNumPoints(OGRGeometryH ogrg) {
	int rings, geos, total;
	int type = wkbFlatten(((OGRGeometry *) ogrg)->getGeometryType());
	switch (type) {
		case wkbPoint:
		case wkbLineString:
			return OGR_G_GetPointCount(ogrg);
		case wkbPolygon:
			total = OGR_G_GetPointCount(((OGRPolygon*)ogrg)->getExteriorRing());
			total += OGR_G_GetPointCount(((OGRPolygon*)ogrg)->getInteriorRing(0));
			return total;
		case wkbMultiPolygon:
		case wkbMultiLineString:
		case wkbMultiPoint:
			geos = ((OGRGeometryCollection*)ogrg)->getNumGeometries();
			total = 0;
			for(int i = 0; i < geos; i++)
				total += OGRGetNumPoints(((OGRGeometryCollection*)ogrg)->getGeometryRef(i));
			return total;
		default:
			printf("Geometry type not count points %d\n", type);
			return 0;
	}
}

void GetEnvelopeCoordinates(Envelope p, double *MinX, double *MinY, double *MaxX, double *MaxY) {
	*MinX = p.MinX;
	*MinY = p.MinY;
	*MaxX = p.MaxX;
	*MaxY = p.MaxY;
}
Envelope GetEnvelope(double MinX, double MinY, double MaxX, double MaxY) {
	return Envelope{MinX, MinY, MaxX, MaxY};
}
