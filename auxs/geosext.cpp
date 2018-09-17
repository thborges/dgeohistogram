/*
 * geosext.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: thborges
 */

#include <cpl_conv.h>
#include <stdarg.h>
#include <assert.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Envelope.h>
#include <geos/io/WKBWriter.h>
#include <geos/io/WKTWriter.h>
#include <geos/io/Writer.h>
#include "geosext.h"

extern "C" const GEOSIEnvelope *GEOSInternalEnvelope(const GEOSGeometry *g1);
extern "C" void GEOSEnvelopeGetXY(const GEOSGeometry *e, double *x1, double *x2, double *y1, double *y2);
extern "C" size_t GEOSGetNumPoints(const GEOSGeometry *e);
extern "C" void geos_messages(const char *fmt, ...);

extern "C" void *GEOSCreateWKBBuffer();
extern "C" void *GEOSCreateWKTBuffer();
extern "C" GEOSWKBWriter *GEOSCreateWKBWritter(GEOSContextHandle_t extHandle);
extern "C" const char *GEOSGeomToWKB(const GEOSWKBWriter *w, const void *buff, const GEOSGeometry *g, int *size);
extern "C" void GEOSDestroyWKBWritter(GEOSWKBWriter *w, void *buff);
extern "C" GEOSWKTWriter *GEOSCreateWKTWritter();
extern "C" const char *DgeoGeomToWKT(const GEOSWKTWriter *w, const void *buff, const GEOSGeometry *g);
extern "C" void GEOSDestroyWKTWritter(GEOSWKTWriter *w, void *buff);

const GEOSIEnvelope *GEOSInternalEnvelope(const GEOSGeometry *g1) {
	const geos::geom::Envelope *e = ((const geos::geom::Geometry*)g1)->getEnvelopeInternal();
    return (GEOSIEnvelope *)e;
}

void GEOSEnvelopeGetXY(const GEOSGeometry *g1, double *x1, double *x2, double *y1, double *y2) {
	const geos::geom::Envelope *e = ((const geos::geom::Geometry*)g1)->getEnvelopeInternal();
	if (e) {
		*x1 = e->getMinX();
		*y1 = e->getMinY();
		*x2 = e->getMaxX();
		*y2 = e->getMaxY();
	}
	else
		*x1 = *x2 = *y1 = *y2 = 0.0;
}

size_t GEOSGetNumPoints(const GEOSGeometry *e) {
	return ((const geos::geom::Geometry*)e)->getNumPoints();
}

void geos_messages(const char *fmt, ...) {
    va_list(ap);
	va_start(ap, fmt);
	fprintf(stderr, "error: GEOS: ");
	vfprintf(stderr, fmt, ap);
	fprintf(stderr, "\n");
	va_end(ap);
}



typedef struct GEOSContextHandleInternal
{
    const geos::geom::GeometryFactory *geomFactory;
    GEOSMessageHandler NOTICE_MESSAGE;
    GEOSMessageHandler ERROR_MESSAGE;
    int WKBOutputDims;
    int WKBByteOrder;
    int initialized;
} GEOSContextHandleInternal_t;

class raw_buffer : public std::streambuf
{
public:
    raw_buffer() {
        buff = (char_type*)malloc(sizeof(char_type)*current_size);
        setp(&buff[0], &buff[current_size]);
    };
    ~raw_buffer() {
    	free(buff);
    }
    int_type overflow(int_type c) override;
    std::streamsize xsputn(const char_type*, std::streamsize) override;
    int sync() override {
        setp(&buff[0], &buff[current_size]);
        return 0;
    };
    const char_type *get_buffer(int *size) { *size = (pptr()-&buff[0])/sizeof(char_type); return buff; };
private:
    int current_size = 16;
    char_type *buff;
};

class raw_ostream : private virtual raw_buffer
                  , public std::ostream
{
public:
    raw_ostream() : raw_buffer(), std::ostream(this) { }

    const char *get_buffer(int *size) {
    	return this->raw_buffer::get_buffer(size);
    }
};

raw_buffer::int_type raw_buffer::overflow(raw_buffer::int_type c) {
    int old_size = current_size;
    current_size += current_size;
    buff = (char_type*)realloc(buff, sizeof(char_type)*current_size);
    if (!buff)
        return traits_type::eof();
    setp(&buff[0], &buff[current_size]);
    pbump(old_size);
    *pptr() = traits_type::to_char_type(c);
    pbump(1);
    return c;
}

std::streamsize raw_buffer::xsputn(const char_type* str, std::streamsize count) {
    for (int i = 0; i < count; ++i) {
        sputc(str[i]);
    }
    return count;
}

void *GEOSCreateWKBBuffer() {
	return new raw_ostream();
}

void *GEOSCreateWKTBuffer() {
	return new std::string();
}

GEOSWKBWriter *GEOSCreateWKBWritter(GEOSContextHandle_t extHandle) {
	GEOSContextHandleInternal_t *handle = 0;
	handle = reinterpret_cast<GEOSContextHandleInternal_t*>(extHandle);
	const geos::io::WKBWriter *w = new geos::io::WKBWriter(2, handle->WKBByteOrder);
	return (GEOSWKBWriter*)w;
}

const char *GEOSGeomToWKB(const GEOSWKBWriter_t *w, const void *buff, const GEOSGeometry *g, int *size) {
	raw_ostream *ostreambuf = ((raw_ostream*)buff);
	ostreambuf->flush();
	((geos::io::WKBWriter*)w)->write(*(const geos::geom::Geometry*)g, *ostreambuf);
	return ostreambuf->get_buffer(size);
}

void GEOSDestroyWKBWritter(GEOSWKBWriter_t *w, void *buff) {
	delete (geos::io::WKBWriter*)w;
	delete (raw_ostream*)buff;
}

GEOSWKTWriter *GEOSCreateWKTWritter() {
	const geos::io::WKTWriter *w = new geos::io::WKTWriter();
	return (GEOSWKTWriter*)w;
}

const char *DgeoGeomToWKT(const GEOSWKTWriter_t *w, const void *buff, const GEOSGeometry *g) {
	/*std::string *sbuff = ((std::string*)buff);
	geos::io::WKTWriter *writerwkt = ((geos::io::WKTWriter*)w);
	const geos::geom::Geometry* geo = (const geos::geom::Geometry*)g;
	*sbuff = writerwkt->write(geo);
	return sbuff->c_str();*/
	return NULL;
}

void GEOSDestroyWKTWritter(GEOSWKTWriter_t *w, void *buff) {
	delete (geos::io::WKTWriter*)w;
	delete (std::string*)buff;
}
