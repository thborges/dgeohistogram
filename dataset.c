/*
 * dataset.c
 *
 *  Created on: 26/07/2014
 *      Author: thborges
 */

#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <ogrext.h>
#include "dataset.h"
#include "uthash.h"
#include "wkbconvert.h"
#include "glibwrap.h"

dataset_page *dataset_getpage(dataset *dh, dataset_segment *seg, unsigned pid);
dataset_page *dataset_alloc_page(unsigned pid);

void dataset_fatal(const char *str) {
	fprintf(stderr, "%s\n", str);
	exit(1);
}

inline __attribute__((always_inline))
dataset_leaf *dataset_get_leaves_addr(dataset *dh, dataset_page *pg, unsigned pos) {
	dataset_leaf *l = &pg->data->leaves[pos * dh->metadata.geocount];
	return l;
}

void dataset_write_page(const dataset *dh, dataset_page *pg) {
	if (pg->wstate == DIRTY) {
		int aux = lseek(dh->fpages, pg->id * DATASET_PAGE_SIZE, SEEK_SET);
		assert(aux != -1);
		aux = write(dh->fpages, pg->memory, DATASET_PAGE_SIZE);
		assert(aux == DATASET_PAGE_SIZE);
		pg->wstate = CLEAN;
		//fprintf(stderr, "Wrote page %s:%d\n", dh->metadata.name, pg->id);
	}
}

dataset_segment *dataset_get_seg(dataset *dh, unsigned sid) {
	dataset_segment *seg = NULL;
	HASH_FIND_INT(dh->segments, &sid, seg);
	return seg;
}

dataset_page *dataset_free_and_reuse_page(dataset *dh, dataset_segment *seg, unsigned newid) {

	// remove the first page == the oldest page
	dataset_page *pg = NULL, *tmp;
	HASH_ITER(hh, seg->page_cache, pg, tmp) {
		if (pg->state == UNLOCKED)
			break;
	}

	if (pg == NULL) { // all pages state == FILLING
		pg = dataset_alloc_page(newid);
		return pg;
	}

	dataset_write_page(dh, pg);
	HASH_DEL(seg->page_cache, pg);

	// setup memory as a new page
	pg->id = newid;
	pg->data->used = 0;
	pg->state = UNLOCKED;
	pg->wstate = DIRTY;

	return pg;
}

dataset_page *dataset_alloc_page(unsigned pid) {
	dataset_page *pg = g_new(dataset_page, 1);
	assert(pg != NULL);

	pg->id = pid;
	pg->state = UNLOCKED;
	pg->wstate = DIRTY;
	pg->memory = g_new0(unsigned char, DATASET_PAGE_SIZE);
	pg->data->used = 0;
	return pg;
}

unsigned dataset_get_new_page_mem(dataset *dh, dataset_segment *seg, dataset_page **pg) G_GNUC_WARN_UNUSED_RESULT;
unsigned dataset_get_new_page_mem(dataset *dh, dataset_segment *seg, dataset_page **pg) {

	pthread_mutex_lock(dh->mem_get_lock);
	unsigned pid = dh->metadata.pagecount;
	dh->metadata.pagecount++;
	pthread_mutex_unlock(dh->mem_get_lock);

	if (DATASET_HASH_FULL(dh, seg))
		*pg = dataset_free_and_reuse_page(dh, seg, pid);
	else
		*pg = dataset_alloc_page(pid);

	(*pg)->state = LOCKED;
	HASH_ADD_INT(seg->page_cache, id, *pg);

	return pid;
}

dataset_page *dataset_get_page_mem(dataset *dh, dataset_segment *seg, unsigned pid) {

	assert(pid < dh->metadata.pagecount);

	dataset_page *pg;
	if (DATASET_HASH_FULL(dh, seg))
		pg = dataset_free_and_reuse_page(dh, seg, pid);
	else
		pg = dataset_alloc_page(pid);

	pg->state = LOCKED;
	HASH_ADD_INT(seg->page_cache, id, pg);

	return pg;
}

void dataset_new_segment_item(dataset *dh, dataset_segment *seg) {
	dataset_segment_item *sitem = g_new0(dataset_segment_item, 1);
	assert(sitem != NULL);

	if (seg->fillingitem) {
		seg->fillingpage->state = UNLOCKED;
		seg->fillingitem->next = sitem;
	}
	sitem->pid = dataset_get_new_page_mem(dh, seg, &seg->fillingpage);
	seg->fillingitem = sitem;
	seg->itemcount++;
}

dataset_segment *dataset_alloc_first_segment(dataset *dh) {
	dataset_segment *seg = g_new0(dataset_segment, 1);
	assert(seg != NULL);

	seg->sid = 0;
	seg->count = 0;
	seg->itemcount = 0;
	dataset_new_segment_item(dh, seg);
	seg->items = seg->fillingitem;

	return seg;
}

dataset_segment *dataset_new_seg(dataset *dh, unsigned sid) {

	dataset_segment *seg = NULL;
	HASH_FIND_INT(dh->segments, &sid, seg);
	if (seg != NULL)
		return seg;

	//printf("New segment %d for dataset %s.\n", sid, dh->metadata.name);

	seg = g_new(dataset_segment, 1);
	assert(seg != NULL);

	seg->sid = sid;
	seg->count = 0;
	seg->itemcount = 0;
	seg->fillingitem = NULL;
	seg->fillingpage = NULL;
	seg->page_cache = NULL;

	dataset_new_segment_item(dh, seg);
	seg->items = seg->fillingitem;

	HASH_ADD_INT(dh->segments, sid, seg);
	HASH_FIND_INT(dh->segments, &sid, seg);
	assert(seg != NULL);

	return seg;
}

void dataset_getpage_disk(dataset *dh, dataset_segment *seg, dataset_page *pg, int fpages) {
	int bytes = 0;
	int repeat = 1000;
    int wait = 10000;

    do {
		int aux = lseek(fpages, pg->id * DATASET_PAGE_SIZE, SEEK_SET);
		assert(aux != -1);
		bytes = read(fpages, pg->memory, DATASET_PAGE_SIZE);
		if (bytes != DATASET_PAGE_SIZE) {
			fprintf(stderr, "Read page %s %d/%d out of bounds. Read bytes %d. Waiting %dms.\n", 
                dh->metadata.name, pg->id, dh->metadata.pagecount, bytes, wait/1000);
			repeat--;
			usleep(wait);
            wait *= 2;
		}
	}
	while (bytes != DATASET_PAGE_SIZE && repeat > 0);

	assert(bytes == DATASET_PAGE_SIZE);
	assert(pg->data->used < DATASET_PAGE_CAPACITY(dh));

	pg->wstate = CLEAN;
}

dataset_page *dataset_getpage(dataset *dh, dataset_segment *seg, unsigned pid) {

	dataset_page *pg;

	HASH_FIND_INT(seg->page_cache, &pid, pg);
	if (pg) {
		dataset_page* tailpage = (dataset_page*)seg->page_cache->hh.tbl->tail->prev;
		if (tailpage && (dataset_page*)tailpage->hh.next != pg) { // move the page to last used
			HASH_DEL(seg->page_cache, pg);
			HASH_ADD_INT(seg->page_cache, id, pg);
		}
		pg->state = LOCKED;
		return pg;
	}

    // in memory datasets all pages are on cache
	assert(dh->metadata.memory_dataset == 0);

	pg = dataset_get_page_mem(dh, seg, pid);
    dataset_getpage_disk(dh, seg, pg, dh->fpages);

	return pg;
}

dataset *dataset_create_full(const char *name, unsigned short geocount, char memory_dataset) {

	dataset *result = g_new0(dataset, 1);
	assert(result != NULL);

	result->mem_get_lock = g_new(pthread_mutex_t, 1);
	int ini = pthread_mutex_init(result->mem_get_lock, NULL);
	assert(ini == 0);

	result->metadata.count = 0;
	result->metadata.geocount = geocount;
	result->metadata.pagecount = 0;
	result->metadata.memory_dataset = memory_dataset;
	result->metadata.hist.mbr = emptymbr;
	result->metadata.hist.xtics = NULL;
	result->metadata.hist.ytics = NULL;
	result->metadata.hist.hcells = NULL;
	result->metadata.name = strdup(name);
	result->metadata.x_average = 0.0;
	result->metadata.y_average = 0.0;
	result->metadata.x_psa = 0.0;
	result->metadata.y_psa = 0.0;
	result->segments = NULL;

	dataset_segment *seg = dataset_alloc_first_segment(result);
	HASH_ADD_INT(result->segments, sid, seg);

	if (!memory_dataset) {
		if (strlen(result->metadata.name) == 0) {
			g_free(result->metadata.name);
			result->metadata.name = g_new(char, DATASET_NAME_MAX);
			sprintf(result->metadata.name, "%lx", (long)result);
		}

		char filename[100];
		sprintf(filename, "./tmp/%s", result->metadata.name);
		result->fpages = open(filename, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
		if (result->fpages == -1)
			perror("Error creating temporary file for dataset");
		assert(result->fpages != -1);
	}
	else
		result->fpages = -1;

	return result;
}

void dataset_sync(dataset *dh) {
	dataset_segment *seg, *seg_tmp;
	HASH_ITER(hh, dh->segments, seg, seg_tmp) {
	dataset_page *page, *page_tmp;
		HASH_ITER(hh, seg->page_cache, page, page_tmp) {
			if (page->wstate != CLEAN)
				dataset_write_page(dh, page);
		}
	}
}

dataset* dataset_clone(dataset *dh) {

	dataset_sync(dh); // force dirty pages write

	dataset *result = g_new(dataset, 1);
	assert(result != NULL);

	*result = *dh;
	result->metadata.cloned = TRUE;
	result->segments = NULL;

	result->mem_get_lock = g_new(pthread_mutex_t, 1);
	int ini = pthread_mutex_init(result->mem_get_lock, NULL);
	assert(ini == 0);
	
	dataset_segment *current, *tmp;
	HASH_ITER(hh, dh->segments, current, tmp) {
		dataset_segment *nseg = g_new0(dataset_segment, 1);
		nseg->sid = current->sid;
		nseg->count = current->count;
		nseg->itemcount = current->itemcount;
		nseg->items = current->items;
		nseg->page_cache = NULL;
		HASH_ADD_INT(result->segments, sid, nseg);
	}

	if (!result->metadata.memory_dataset) {
		char filename[100];
		sprintf(filename, "./tmp/%s", dh->metadata.name);
		result->fpages = open(filename, O_RDONLY);
		if (result->fpages == -1)
			perror("Error reopening temporary file to clone dataset");
	}

	return result;
}

void dataset_seg_sync(dataset *dh, dataset_segment *seg) {
	// sync the pages of the seg on the original dataset
    pthread_mutex_lock(dh->mem_get_lock);

	dataset_segment_item *it = seg->items;
	while (it != NULL) {
		dataset_page *pg = NULL;

		HASH_FIND_INT(seg->page_cache, &it->pid, pg);
		if (pg)
			dataset_write_page(dh, pg);

		it = it->next;
	}

	pthread_mutex_unlock(dh->mem_get_lock);
}

dataset* dataset_clone_seg(dataset *dh, dataset_segment *seg, dataset_segment **clone_seg) {

	dataset *result = g_new(dataset, 1);
	assert(result != NULL);

	*result = *dh;
	result->metadata.cloned = TRUE;
	result->segments = NULL;

	result->mem_get_lock = g_new(pthread_mutex_t, 1);
	int ini = pthread_mutex_init(result->mem_get_lock, NULL);
	assert(ini == 0);

	dataset_segment *current = seg;
	assert(current != NULL);

	// add the sid specified
	dataset_segment *nseg = g_new0(dataset_segment, 1);
	nseg->sid = current->sid;
	nseg->count = current->count;
	nseg->itemcount = current->itemcount;
	nseg->items = current->items;
	nseg->page_cache = NULL;
	HASH_ADD_INT(result->segments, sid, nseg);
	*clone_seg = nseg;

	// sync the pages of the seg on the original dataset
	dataset_seg_sync(dh, seg);

	if (!result->metadata.memory_dataset) {
		char filename[100];
		sprintf(filename, "./tmp/%s", dh->metadata.name);
		result->fpages = open(filename, O_RDONLY);
		if (result->fpages == -1)
			perror("Error reopening temporary file to clone dataset");
	}

	return result;
}

dataset* dataset_create(const char *name, unsigned short geocount) {
	return dataset_create_full(name, geocount, FALSE);
}

long long dataset_destroy_seg_full(dataset *dh, dataset_segment *cseg, char destroy_geos) {
    long long memused = 0;
    
    // release geo objects
    if (destroy_geos && !dh->metadata.cloned) {
        dataset_iter di;
       	dataset_foreach_seg(di, dh, cseg) {
            for (int i = 0; i < dh->metadata.geocount; i++) {
	         	GEOSPreparedGeom_destroy(di.item[i].pgeo);
		        GEOSGeom_destroy(di.item[i].geo);
            }
   	    }
    }

	// release page cache of the segment
	dataset_page *page, *page_tmp;
	HASH_ITER(hh, cseg->page_cache, page, page_tmp) {
		HASH_DEL(cseg->page_cache, page);
		g_free(page->memory);
		g_free(page);
        memused += sizeof(dataset_page) + DATASET_PAGE_SIZE;
    }

    // release segment items
    if (!dh->metadata.cloned) {
		dataset_segment_item *i = cseg->items;
   		dataset_segment_item *next;
    	do {
	        next = i->next;
		    g_free(i);
			i = next;
            memused += sizeof(dataset_segment_item);
        } while (i != NULL);
    }

    // release segment and hash
	HASH_DEL(dh->segments, cseg);
	g_free(cseg);

    return memused;
}

void dataset_destroy_seg(dataset *dh, dataset_segment *seg) {
    dataset_destroy_seg_full(dh, seg, TRUE);
}

void dataset_destroy_full(dataset *dh, char destroy_geos) {
	if (!dh)
		return;

	long long memused = 0;

	dataset_segment *cseg, *tmpseg;
	HASH_ITER(hh, dh->segments, cseg, tmpseg) {
        memused += dataset_destroy_seg_full(dh, cseg, destroy_geos);
        memused += sizeof(dataset_segment);
	}

    // release dataset and histogram stuff
	if (!dh->metadata.cloned) {
        fprintf(stderr, "Dataset %s mem usage: ", dh->metadata.name);

		g_free(dh->metadata.name);
        memused += DATASET_NAME_MAX;

		// histogram stuff
		g_free(dh->metadata.hist.xtics);
		g_free(dh->metadata.hist.ytics);
		g_free(dh->metadata.hist.hcells);
        memused += 
            sizeof(double) * (dh->metadata.hist.xsize + dh->metadata.hist.ysize) +
            sizeof(histogram_cell) * dh->metadata.hist.xsize * dh->metadata.hist.ysize;
	}
		
    if (!dh->metadata.memory_dataset)
        close(dh->fpages);

	pthread_mutex_destroy(dh->mem_get_lock);
    g_free(dh->mem_get_lock);
    
    memused += sizeof(pthread_mutex_t);
    memused += sizeof(dataset);
    if (!dh->metadata.cloned)
        fprintf(stderr, "%lli\n", memused);

	g_free(dh);
}

void dataset_destroy(dataset *dh) {
	dataset_destroy_full(dh, TRUE);
}

dataset* dataset_create_mem(const char *name, unsigned short geocount) {
	return dataset_create_full(name, 1, TRUE);
}

dataset_leaf *dataset_add_seg(dataset *dh, dataset_segment *seg) {

	if (seg == NULL) printf("Segment %s_%d_%d is null\n", dh->metadata.name, (unsigned short)seg->sid, seg->sid<<16);
	assert(seg != NULL);

	unsigned next = (seg->fillingpage->data->used+1);
	if (next >= DATASET_PAGE_CAPACITY(dh)) { // the space on the page is gone
		dataset_new_segment_item(dh, seg);
	}

	dataset_leaf *l = dataset_get_leaves_addr(dh, seg->fillingpage, seg->fillingpage->data->used);
	seg->fillingpage->data->used++;
	seg->count++;
	dh->metadata.count++;

	return l;

}

dataset_leaf *dataset_add(dataset *dh) {
	return dataset_add_seg(dh, dataset_get_seg(dh, 0));
}

void dataset_concat(dataset *dh, dataset *di) {
/*	if (dh->metadata.geocount != di->metadata.geocount)
		fprintf(stderr, "Dataset's to concatenate are of different size.\n");

	dataset_segment *seg = dataset_get_seg(dh, 0);
	assert(seg != NULL);

	// move pages from di
	dataset_page *dipg;
	for(unsigned i = 0; i < di->metadata.pagecount; i++) {
		dipg = dataset_getpage(di, i);
		HASH_DEL(di->pages, dipg);

		dipg->id = dh->metadata.pagecount;
		dh->metadata.pagecount++;

		dipg->wstate = DIRTY;
		HASH_ADD_INT(dh->pages, id, dipg);

		dataset_segment_item *sitem = g_new(dataset_segment_item, 1);
		assert(sitem != NULL);

		sitem->pid = dipg->id;
		sitem->next = NULL;

		seg->fillingitem->next = sitem;
		seg->fillingitem = sitem;
		seg->itemcount++;

		dh->metadata.count += dipg->data->used;
		seg->count += dipg->data->used;
	}
	seg->fillingpage = dipg;
	seg->fillingpage->state = LOCKED;

	di->metadata.pagecount = 0;
	dataset_destroy_full(di, FALSE);
*/
}

dataset_iter_seg dataset_first_seg(dataset *dh, dataset_segment *seg) {
	assert(seg != NULL);

	dataset_iter_seg i;
	i.dh = dh;
	i.seg = seg;
	i.si = seg->items;
	i.position = 0;
    i.fpages = -1;
    i.page = NULL;
    i.item = NULL;

    if (dh->metadata.memory_dataset)
        i.page = dataset_getpage(dh, seg, i.si->pid);
    else {
        // each iterator open the file because, in general,
        // a iterator is used inside a thread on distributed dgeo
        char filename[7+strlen(dh->metadata.name)];
		sprintf(filename, "./tmp/%s", dh->metadata.name);
		i.fpages = open(filename, O_RDONLY);
		if (i.fpages == -1) {
			perror("Error reopening temporary file to iter on dataset");
            return i;
        }

        // sync the pages because the iterator always get the pages from disk
        dataset_seg_sync(dh, seg);

    	i.page = dataset_alloc_page(i.si->pid);
        dataset_getpage_disk(dh, seg, i.page, i.fpages);
    }

    dataset_fetch_seg(&i);
	return i;
}

inline __attribute__((always_inline))
void dataset_release_iter(dataset_iter_seg *i) {
    if (i->page) {
    	i->page->state = UNLOCKED;
        if (!i->dh->metadata.memory_dataset) {
            g_free(i->page->memory);
            g_free(i->page);
            close(i->fpages);
        }
    	i->page = NULL;
	    i->item = NULL;
    }
}

inline __attribute__((always_inline))
void dataset_fetch_seg(dataset_iter_seg *i) {
	if (i->page->data->used > i->position)
		i->item = dataset_get_leaves_addr(i->dh, i->page, i->position);
	else {
		// this condition only occurs when segment was started, but has no item
		// otherwise, dataset_next_seg will return before call this procudure.
		assert(i->position == 0);
		dataset_release_iter(i);
	}
}


void dataset_next_seg(dataset_iter_seg *i) {
	unsigned next = (i->position+1);
	if (next == i->page->data->used) { // the last item on the page
		i->page->state = UNLOCKED;
		dataset_segment_item *si = i->si->next;
		if (si) {
			i->si = si;
            if (!i->dh->metadata.memory_dataset) {
			    i->page->id = i->si->pid;
                dataset_getpage_disk(i->dh, i->seg, i->page, i->fpages);
            }
            else
                i->page = dataset_getpage(i->dh, i->seg, i->si->pid);
			i->position = 0;
		}
		else {
			dataset_release_iter(i);
			return;
		}
	}
	else // has more items on the page
		i->position++;

	dataset_fetch_seg(i);
}

void dataset_fill_leaf_ext(dataset_leaf *leaf, int index, long long gid, bool cloned, GEOSGeometryH ggeo, Envelope *e, GEOSPreparedGeometry *pgeo) {
	dataset_fill_leaf(leaf, index, ggeo, gid, e, pgeo);
	leaf[index].cloned = cloned;
}

void dataset_fill_leaf_id(dataset_leaf *leaf, int index, long long gid, Envelope *e) {
	dataset_fill_leaf(leaf, index, NULL, gid, e, NULL);
}

void dataset_fill_leaf(dataset_leaf *leaf, int index, GEOSGeometryH ggeo, long long gid, Envelope *e, GEOSPreparedGeometry *pgeo) {
	leaf[index].cloned = false;
	leaf[index].gid = gid;
	leaf[index].geo = ggeo;
	leaf[index].pgeo = pgeo;
	if (e) leaf[index].mbr = *e;
	else GEOSGetEnvelope(ggeo, &leaf[index].mbr);
}

void dataset_set_histogram(dataset *dh, dataset_histogram *hist) {
	dh->metadata.hist = *hist;
	// currently, these are invalid references from another remote histogram
	dh->metadata.hist.xtics = NULL;
	dh->metadata.hist.ytics = NULL;
	dh->metadata.hist.hcells = NULL;
}

dataset_histogram *dataset_get_histogram(dataset *ds) {
	return &ds->metadata.hist;
}

GEOSGeometryH dataset_get_leaf_geo(dataset *dh, dataset_leaf *leaf) {
	if (leaf->gid == -1)
		return leaf->geo;

	assert(dh->temp_ogr_layer && "Missing temporary layer to load geometry.");

	OGRFeatureH f = OGR_L_GetFeature(dh->temp_ogr_layer, leaf->gid);
	OGRGeometryH g = OGR_F_GetGeometryRef(f);
	GEOSGeometryH geosg = convertOGRToGEOS(g);
	OGR_F_Destroy(f);
	return geosg;
}

