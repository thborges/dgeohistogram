/*
 * wkbconvert.c
 *
 *  Created on: Jul 9, 2015
 *      Author: thborges
 */

#include <wkbconvert.h>
#include <assert.h>

size_t _wkb_current_size = 0;
unsigned char *_wkb_buffer = NULL;

GEOSGeometry *convertOGRToGEOS(OGRGeometryH *geo) {
	size_t wkb_size = OGR_G_WkbSize(geo);
	if (wkb_size > _wkb_current_size){
		unsigned char* pNewBuffer = (unsigned char*)realloc(_wkb_buffer, sizeof(unsigned char) * wkb_size);
        if(pNewBuffer){
            _wkb_buffer = pNewBuffer;
            _wkb_current_size = wkb_size;
        }else{
            if(_wkb_buffer){
                free(_wkb_buffer);
                _wkb_buffer = NULL;
            }
            _wkb_current_size = 0;
            assert(0 && "Out of Memory or Memory Fragmentation Issue !");
        }
    }

	OGR_G_ExportToWkb(geo, (OGRwkbByteOrder)(htonl(1) == 1 ? 0 : 1), _wkb_buffer);
	GEOSGeometry *ggeo = GEOSGeomFromWKB_buf(_wkb_buffer, wkb_size);
	return ggeo;
}
