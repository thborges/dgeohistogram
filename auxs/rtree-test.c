#include <assert.h>
#include "rtree-test.h"
#include "ogrext.h"

void count_levels(rtree_node *node, int *levels) {
    if (node->type == DIRECTORY) {
    	count_levels(node->dirs[0], levels);
    }
    *levels += 1;
}

count_node_struct count_node(rtree_node *node, const int m) {
	count_node_struct c;
    c.leaves = c.dirs = c.items = c.subleaves = 0;
    
    if (node->type == DIRECTORY) {
		for(int i=0; i < node->used; i++) {
			count_node_struct caux = count_node(node->dirs[i], m);
            c.items += caux.items;
            c.dirs += caux.dirs;
            c.leaves += caux.leaves;
            c.subleaves += caux.subleaves;
        }
		c.dirs += 1;
	}
	else {
		if (node->used < minm(m))
			c.subleaves += .5;
		else
			c.leaves += 1;
        c.items += node->used;
	}
	return c;
}

count_node_struct count_nodes(rtree_root *rtree) {
	count_node_struct c = count_node(rtree->root, rtree->m);
	return c;
}

void check_mbrs_recursive(rtree_node *node) {
	if (node->type == DIRECTORY)
		for(int i = 0; i < node->used; i++)
			check_mbrs_recursive(node->dirs[i]);
#ifndef NDEBUG
	Envelope e = node->mbr;
 	Envelope f = rtree_compute_mbr(node);
#endif
	if (!(e.MinX == f.MinX && e.MinY == f.MinY && e.MaxX == f.MaxX && e.MaxY == f.MaxY)) {
		printf("%f %f, %f %f\n", e.MinX, e.MinY, e.MaxX, e.MaxY);
		printf("%f %f, %f %f\n", f.MinX, f.MinY, f.MaxX, f.MaxY);
	}
	assert(e.MinX == f.MinX && e.MinY == f.MinY && e.MaxX == f.MaxX && e.MaxY == f.MaxY);
}

void check_mbrs(rtree_root *rtree) {
	check_mbrs_recursive(rtree->root);
}

OGRGeometryH geometryFromEnvelope(const Envelope env, OGRSpatialReferenceH refspatial) {
	char wkt[512];
	sprintf(wkt, "POLYGON((%e %e, %e %e, %e %e, %e %e, %e %e))",
		env.MinX, env.MinY,
		env.MaxX, env.MinY,
		env.MaxX, env.MaxY,
		env.MinX, env.MaxY,
		env.MinX, env.MinY);
	OGRGeometryH geo;
	char *wktp = &wkt[0];
	OGR_G_CreateFromWkt(&wktp, refspatial, &geo);
	return geo;
}

void calc_area_recursive(rtree_node *node, double *areaD, GList **lmbrs) {
	if (node->type == DIRECTORY) {
		for(int i = 0; i < node->used; i++)
			calc_area_recursive(node->dirs[i], areaD, lmbrs);
	
		if (node->dirs[0]->type == LEAF) {
			for(int i = 0; i < node->used; i++) {
				*lmbrs = g_list_prepend(*lmbrs, &node->dirs[i]->mbr);
			}
			*areaD += ENVELOPE_AREA(node->mbr);
		}
		else
		for(int i = 0; i < node->used; i++) {
			*areaD += ENVELOPE_AREA(node->dirs[i]->mbr);
		}
	}
}

void print_overlap(rtree_root *rtree, OGRSpatialReferenceH refspatial) {
	setlocale(LC_NUMERIC, "C");

	double areaD = 0.0;
	int levels = 0;
	count_levels(rtree->root, &levels);

	GList *lmbrs = NULL;
	calc_area_recursive(rtree->root, &areaD, &lmbrs);
	double leaves_overlap = 0.0;
	for (GList *ni = g_list_next(lmbrs);
		ni != NULL;
		ni = g_list_next(ni)) {
		Envelope *a = (Envelope *)ni->data;
		for(GList *cp = g_list_next(ni); cp != NULL; cp = g_list_next(cp)) {
			Envelope *b = (Envelope *)cp->data;
			leaves_overlap += ENVELOPE_AREA(EnvelopeIntersection(*a, *b));
		}
	}

	printf("%10.4f|%10.4f|", (areaD / (levels-1)) / ENVELOPE_AREA(rtree->root->mbr), leaves_overlap);
}

void fill_layers_for_levels(const rtree_node *node, const rtree_node *parent, const int level, const OGRLayerH *layers) {
	int i;	

	if (node->type == DIRECTORY) {
		for(i = 0; i < node->used; i++) {
			fill_layers_for_levels(node->dirs[i], node, level+1, layers);
		}
	}

	for(i = 0; i < node->used; i++) {
		int usedentries;
		OGRGeometryH geo = NULL;
		if (node->type == LEAF) {
			unsigned char *wkb_buff;
			size_t size;
			wkb_buff = GEOSGeomToWKB_buf(node->leaves[i].geo, &size);
			OGR_G_CreateFromWkb(wkb_buff, NULL, &geo, (int)size);
			usedentries = 0;
		}
		else {
			geo = geometryFromEnvelope(node->dirs[i]->mbr, NULL);
			usedentries = node->dirs[i]->used;
		}

		OGRFeatureH feature = OGR_F_Create(OGR_L_GetLayerDefn(layers[level]));
		OGR_F_SetFieldInteger(feature, OGR_F_GetFieldIndex(feature, "ParentID"), parent != NULL? parent->index : -1);
		OGR_F_SetFieldInteger(feature, OGR_F_GetFieldIndex(feature, "Used"), usedentries);
		(node->type == DIRECTORY) ? OGR_F_SetGeometryDirectly(feature, geo) : OGR_F_SetGeometry(feature, geo);
		OGR_L_CreateFeature(layers[level], feature);
		OGR_F_Destroy(feature);

		if (node->type == LEAF)
			OGR_G_DestroyGeometry(geo);
	}

}

void write_dirs_shapefile(rtree_root *rtree) {
	const char *driverName = "ESRI Shapefile";
	OGRSFDriverH driver;
	
	setlocale(LC_NUMERIC, "C");

	driver = OGRGetDriverByName(driverName);
	if (!driver) {
		fprintf(stderr, "Shapefile driver not found");
		return;
	}

	OGR_Dr_DeleteDataSource(driver, "Level.0.shp");
	OGRDataSourceH ds = OGR_Dr_CreateDataSource(driver, "Level.0.shp", NULL);
	if (!ds) {
		fprintf(stderr, "Cannot create boundbox.shp file.");
		return;
	}


	// type of leaf geom
	rtree_node *n = rtree->root;
	while (n->type != LEAF)
		n = n->dirs[0];
	int geotype = GEOSGeomTypeId(n->leaves[0].geo);
	OGRwkbGeometryType ogrtypeleaf = wkbMultiPolygon;
	switch (geotype) {
		case GEOS_POINT: 
			ogrtypeleaf = wkbPoint; 
			break;
		case GEOS_LINESTRING:
			ogrtypeleaf = wkbLineString;
			break;
	}

	int height = rtree_height(rtree);
	OGRFieldDefnH fielddef, fielddef2;
	OGRLayerH *layers = g_new(OGRLayerH, height);
	for(int l=0; l<height; l++) {
		char s[20];
		sprintf(s, "Level.%d", l);

		if (l == height-1) // leaf level
			layers[l] = OGR_DS_CreateLayer(ds, s, NULL, ogrtypeleaf, NULL);
		else
			layers[l] = OGR_DS_CreateLayer(ds, s, NULL, wkbMultiPolygon, NULL);


		if (!layers[l]) {
			fprintf(stderr, "Cannot create layer for level %d.", l);
			return;
		}

		fielddef = OGR_Fld_Create("ParentID", OFTInteger);
		if (!fielddef || (OGRERR_NONE != OGR_L_CreateField(layers[l], fielddef, TRUE))) {
			fprintf(stderr, "Cannot create field ParentID for level %d.", l);
			return;
		}
		OGR_Fld_Destroy(fielddef);

		fielddef2 = OGR_Fld_Create("Used", OFTInteger);
		if (!fielddef2 || (OGRERR_NONE != OGR_L_CreateField(layers[l], fielddef2, TRUE))) {
			fprintf(stderr, "Cannot create field Used for level %d.", l);
			return;
		}
		OGR_Fld_Destroy(fielddef2);
	}


	fill_layers_for_levels(rtree->root, NULL, 0, layers);

	g_free(layers);
	OGR_DS_Destroy(ds);
	
	setlocale(LC_NUMERIC, "");
}


