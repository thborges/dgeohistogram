/*
 * wkbconvert.h
 *
 *  Created on: Jul 9, 2015
 *      Author: thborges
 */

#ifndef UTILS_WKBCONVERT_H_
#define UTILS_WKBCONVERT_H_

#include <geos_c.h>
#include <ogr_api.h>
#include <arpa/inet.h>

GEOSGeometry* convertOGRToGEOS(OGRGeometryH *geo);

#endif /* UTILS_WKBCONVERT_H_ */
