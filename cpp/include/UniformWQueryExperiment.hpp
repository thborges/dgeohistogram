/*
 * Random window query experiment
 *
 *  Created on: 2022-02-16
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include "GenericExperiment.hpp"
#include "SpatialHistogram.hpp"
#include "Dataset.hpp"

#include "rtree.h"
#include "rtree-star.h"
#include "rtree-lazy.h"

// global, to be used to debug precision errors
//#define DEBUG_QUERYNO
#ifdef DEBUG_QUERYNO
int global_qryno;
double global_query_size;
int global_debug_qryno = 74;
double global_debug_qrysize = 0.1;
rtree_root *global_search_rtree;
#endif

class UniformWQueryExperiment: public GenericExperiment {
public:
    UniformWQueryExperiment(Dataset& d, std::vector<SpatialHistogram*>& hs, 
		std::vector<double> qsizes): ds(d), hists(hs), 
		query_sizes(qsizes) {}

    const virtual void run() {
		// cria uma r*
		rtree_root *rtree = NULL;
		rtree = rtree_new_r0(30, 10);
		
		#ifdef DEBUG_QUERYNO
		global_search_rtree = rtree;
		#endif

		std::printf("Building R0-Tree\n");
		for(const DatasetEntry& de: ds.geoms()) {
			EnvelopeC env{de.mbr.MinX, de.mbr.MinY, de.mbr.MaxX, de.mbr.MaxY};
			rtree_append(rtree, de.geo, env);
		}

		std::printf("\nSize\tARE\t\tSTD\t\tSUM\tMethod    \tName\tWorst qry\tError\n");
		const DatasetMetadata& meta = ds.metadata();
		rtree_window_stat stats;
		double width = meta.mbr.MaxX - meta.mbr.MinX;
		double height = meta.mbr.MaxY - meta.mbr.MinY;

		for(double query_size: query_sizes) {
			#ifdef DEBUG_QUERYNO
			global_query_size = query_size;
			#endif

			// metrics for each histogram
			int hs = hists.size();
			double sum_ei[hs];
			double sum_ri[hs];
			double mean[hs];
			double M2[hs];
			double sum_error[hs];
			int worst_queryno[hs];
			int worst_queryerror[hs];
			for(int i = 0; i < hs; i++) {
				sum_ei[i] = 0.0;
				sum_ri[i] = 0.0;
				mean[i] = 0.0;
				M2[i] = 0.0;
				sum_error[i] = 0.0;
				worst_queryno[i] = -1;
				worst_queryerror[i] = INT_MIN;
			}

			double wsize = width * query_size;
			double hsize = height * query_size;

			// quantidade de consultas para cobrir o dataset, considerando
			// uma distribuicao uniforme
			int qtd = ceil(width / wsize) * ceil(height/hsize);
			int qtdqx = ceil(width / wsize);

			//std::printf("Query count: %d, w %f, h %f, w_size %f, h_size %f\n", qtd, width, height, wsize, hsize);

			int n = 0;
			int qryno = 0;
			
			#ifdef DEBUG_QUERYNO
			global_qryno = 0;
			#endif

			while (n < qtd) {
				n++;

				Envelope query;
				query.MinX = meta.mbr.MinX + (qryno / qtdqx) * wsize;
				query.MinY = meta.mbr.MinY + (qryno % qtdqx) * hsize;
				query.MaxX = query.MinX + wsize;
				query.MaxY = query.MinY + hsize;
				EnvelopeC queryc = query;
				
				GEOSGeometryH geoquery = query.toGEOSGeom();

				memset(&stats, 0, sizeof(rtree_window_stat));
				GList *results = rtree_window_search(rtree, geoquery, queryc, &stats);

				// real cardinality from rtree
				int riq = g_list_length(results);

				for(int i = 0; i < hs; i++) {
					// histogram estimate cardinality
					int rhq = hists[i]->estimateWQuery(query);
					
					#ifdef DEBUG_QUERYNO
					if (global_qryno == global_debug_qryno && global_query_size == global_debug_qrysize) {
						char wkt[512];
						std::sprintf(wkt, "POLYGON((%lf %lf, %lf %lf, %lf %lf, %lf %lf, %lf %lf))",
							query.MinX, query.MinY,
							query.MaxX, query.MinY,
							query.MaxX, query.MaxY,
							query.MinX, query.MaxY,
							query.MinX, query.MinY);
						printf("%s\n", wkt);
						printf("%d %-10s\treal: %5d\test: %5d\terr: %5d\n", qryno, 
							hists[i]->name().c_str(), riq, rhq, rhq - riq);
					}
					#endif

					int error = abs(rhq-riq);
					if (error > worst_queryerror[i]) {
						worst_queryerror[i] = error;
						worst_queryno[i] = qryno;
					}

					// average relative error
					sum_ei[i] += error;
					sum_ri[i] += riq;

					// stdev of error
					double delta = error - mean[i];
					mean[i] += delta/(double)n;
					M2[i] += delta*(error - mean[i]);

					// sum error
					sum_error[i] += error;
				}

				g_list_free(results);
				GEOSGeom_destroy(geoquery);
				qryno++;
				
				#ifdef DEBUG_QUERYNO
				global_qryno++;
				#endif
			}
		
			for(int i = 0; i < hs; i++) {
				std::printf("%3.2lf\t%lf\t%lf\t%i\t%-10s\t%s\t%d\t%d\n",
					query_size, 
					sum_ei[i] / (double)sum_ri[i],
					sqrt(M2[i]/(double)n),
					(int)sum_error[i],
					hists[i]->name().c_str(),
					ds.metadata().name.c_str(),
					worst_queryno[i],
					worst_queryerror[i]);
			}
			std::printf("\n");
		}
	}

private:
	Dataset& ds;
    std::vector<SpatialHistogram*>& hists;
	std::vector<double> query_sizes;
};
