/*
 * wkbconvert.c
 *
 *  Created on: Jul 9, 2015
 *      Author: thborges
 */

#include <wkbconvert.h>

size_t _wkb_current_size = 0;
unsigned char *_wkb_buffer = NULL;

GEOSGeometry *convertOGRToGEOS(OGRGeometryH *geo) {
	size_t wkb_size = OGR_G_WkbSize(geo);
	if (wkb_size > _wkb_current_size)
		_wkb_buffer = realloc(_wkb_buffer, sizeof(unsigned char) * wkb_size);

	OGR_G_ExportToWkb(geo, (OGRwkbByteOrder)(htonl(1) == 1 ? 0 : 1), _wkb_buffer);
	GEOSGeometry *ggeo = GEOSGeomFromWKB_buf(_wkb_buffer, wkb_size);
	return ggeo;
}
