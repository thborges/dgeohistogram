#include "ogrext.h"

//EnvelopeC emptymbr = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX};
EnvelopeC emptymbr = {0, 0, 0, 0};

const EnvelopeC EnvelopeIntersection(const EnvelopeC r, const EnvelopeC s) {
    if (ENVELOPE_INTERSECTS(r, s)) {
        EnvelopeC result;
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
const EnvelopeC EnvelopeIntersection2(const EnvelopeC r, const EnvelopeC s) {
    EnvelopeC result = {
    		r.MinX > s.MinX ? r.MinX : s.MinX,
    		r.MinY > s.MinY ? r.MinY : s.MinY,
    		r.MaxX < s.MaxX ? r.MaxX : s.MaxX,
    		r.MaxY < s.MaxY ? r.MaxY : s.MaxY};
    return result;
}

void print_geojson_mbr_local(const EnvelopeC e, const char *id) {
	fprintf(stderr, "{'type': 'Feature', 'geometry': {'type': 'Polygon', 'coordinates': [[");
	fprintf(stderr, "[%f, %f],", e.MinX, e.MinY);
	fprintf(stderr, "[%f, %f],", e.MaxX, e.MinY);
	fprintf(stderr, "[%f, %f],", e.MaxX, e.MaxY);
	fprintf(stderr, "[%f, %f],", e.MinX, e.MaxY);
	fprintf(stderr, "[%f, %f]",  e.MinX, e.MinY);
	fprintf(stderr, "]]}, 'properties': {'name': '%s'}},\n", id);
}

EnvelopeC OGRGetEnvelope(OGRGeometryH ogrg) {
	OGREnvelope ev;
	OGR_G_GetEnvelope(ogrg, &ev);
	EnvelopeC result;
	result.MinX = ev.MinX;
	result.MinY = ev.MinY;
	result.MaxX = ev.MaxX;
	result.MaxY = ev.MaxY;
	return result;
}

void GetEnvelopeCoordinates(EnvelopeC p, double *MinX, double *MinY, double *MaxX, double *MaxY) {
	*MinX = p.MinX;
	*MinY = p.MinY;
	*MaxX = p.MaxX;
	*MaxY = p.MaxY;
}
EnvelopeC GetEnvelope(double MinX, double MinY, double MaxX, double MaxY) {
	EnvelopeC result;
	result.MinX = MinX;
	result.MinY = MinY;
	result.MaxX = MaxX;
	result.MaxY = MaxY;
	return result;
}

void geos_messages(const char *fmt, ...) {
    va_list(ap);
	va_start(ap, fmt);
	fprintf(stderr, "error: GEOS: ");
	vfprintf(stderr, fmt, ap);
	fprintf(stderr, "\n");
	va_end(ap);
}
