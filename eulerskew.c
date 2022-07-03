
//#include "eulerskew.h"
// possibilidade: mover todas as funções do histograma de euler e minskew, os structs do eulerskew é a união dos structs de euler e minskew
//. em fila
// eulerskew_face substituido por eulerskew_face, ambos possuem as mesmas variaveis que são utilizadas pela função
#include <float.h>
#include "histogram.h"
#include "eulerskew.h"
// Utilizar as funções do histograma minskew pra gerar o histograma e particionar os buckets como desejado, depois do histograma minskew ter gerado, utilizar a função
// foreach pra percorrer todos os buckets de minskewlist, em cada bucket, será atribuido valores para os pontos, arestas e faces do histograma de
// euler, de acordo com os pontos de maxX minx maxY e minY do envelope do bucket,  e assim fazer a contagem das intersecções dos objetos com o bucket.
void eulerskew_get_ini_fim(dataset_histogram *dh, Envelope ev, int *xini, int *xfim, int *yini, int *yfim)
{
  *xini = (ev.MinX - dh->mbr.MinX) / dh->xsize;
  *xfim = (ev.MaxX - dh->mbr.MinX) / dh->xsize;
  *yini = (ev.MinY - dh->mbr.MinY) / dh->ysize;
  *yfim = (ev.MaxY - dh->mbr.MinY) / dh->ysize;
  const double epsilon = 1e-100;
  if (ev.MaxX - dh->xtics[*xfim]<epsilon && * xfim> 0)
  {
    (*xfim)--;
  }
  if (ev.MaxY - dh->ytics[*yfim]<epsilon && * yfim> 0)
  {
    (*yfim)--;
  }
  if (*xfim < *xini)
    *xini = *xfim;
  if (*yfim < *yini)
    *yini = *yfim;
}
variance_result_eulerskew eulerskew_calculate_skew_row_col(dataset_histogram *dh, int xini, int xfim, int yini, int yfim)
{
  variance_result_eulerskew r;
  r.n = 0;
  r.mean = 0.0;
  double M2 = 0.0;
  for (int i = xini; i <= xfim; i++)
  {
    for (int j = yini; j <= yfim; j++)
    {
      r.n++;
      double v = GET_HISTOGRAM_CELL(dh, i, j)->cardin;
      double delta = v - r.mean;
      r.mean += delta / (double)r.n;
      M2 += delta * (v - r.mean);
    }
  }
  r.variance = r.n < 2 ? 0.0 : M2 / (double)r.n;
  return r;
}
void eulerskew_calculate_bucket_with_mbr(dataset_histogram *dh, eulerskew_face *b)
{
  int xini, xfim, yini, yfim;
  eulerskew_get_ini_fim(dh, b->mbr, &xini, &xfim, &yini, &yfim);
  variance_result_eulerskew r = eulerskew_calculate_skew_row_col(dh, xini, xfim, yini, yfim);
  b->skew = r.n * r.variance;
  b->cardin = r.n * r.mean;
  b->skew_reduction = NAN;
}
void eulerskew_calculate_skew_reduction(dataset_histogram *dh, eulerskew_face *bucket)
{
  int xini, xfim, yini, yfim;
  eulerskew_get_ini_fim(dh, bucket->mbr, &xini, &xfim, &yini, &yfim);
  // int xqtd = xfim-xini+1;
  // int yqtd = yfim-yini+1;
  bucket->skew_reduction = 0;
  eulerskew_face aux_bucket1, aux_bucket2;
  for (int x = xini; x < xfim; x++)
  {
    // divide on x
    aux_bucket1.mbr = bucket->mbr;
    aux_bucket2.mbr = bucket->mbr;
    aux_bucket1.mbr.MaxX = aux_bucket2.mbr.MinX = dh->xtics[x + 1];
    eulerskew_calculate_bucket_with_mbr(dh, &aux_bucket1);
    eulerskew_calculate_bucket_with_mbr(dh, &aux_bucket2);
    double new_skew = aux_bucket1.skew + aux_bucket2.skew;
    double reduction = bucket->skew - new_skew;
    if (bucket->skew_reduction < reduction)
    {
      bucket->skew_reduction = reduction;
      bucket->split_axis = 'x';
      bucket->split_point = x + 1;
    }
  }
  for (int y = yini; y < yfim; y++)
  {
    // divide on y
    aux_bucket1.mbr = bucket->mbr;
    aux_bucket2.mbr = bucket->mbr;
    aux_bucket1.mbr.MaxY = aux_bucket2.mbr.MinY = dh->ytics[y + 1];
    eulerskew_calculate_bucket_with_mbr(dh, &aux_bucket1);
    eulerskew_calculate_bucket_with_mbr(dh, &aux_bucket2);
    double new_skew = aux_bucket1.skew + aux_bucket2.skew;
    double reduction = bucket->skew - new_skew;
    if (bucket->skew_reduction < reduction)
    {
      bucket->skew_reduction = reduction;
      bucket->split_axis = 'y';
      bucket->split_point = y + 1;
    }
  }
}
minskewLists *eulerskew_generate_hist(dataset *ds, int buckets_num)
{
  minskewLists *listaEulerskew = g_new0(minskewLists, 1);
  dataset_histogram *dh = &ds->metadata.hist;
  eulerskew_face *first_bucket = g_new0(eulerskew_face, 1);
  /*eulerskew_edge  *first_Edge = g_new0(eulerskew_edge, 1);
  eulerskew_vertex *first_Vertex = g_new0(eulerskew_vertex, 1);*/
  first_bucket->mbr = ds->metadata.hist.mbr; // o primeiro bucket será as extremidades do histograma minskew
  /*first_Edge->mbr.MinX = ds->metadata.hist.mbr.MinX;
    first_Edge->mbr.MinY = ds->metadata.hist.mbr.MinY;
    first_Edge->mbr.MaxX = ds->metadata.hist.mbr.MaxX;
    first_Edge->mbr.MaxY = ds->metadata.hist.mbr.MaxY;
    first_Vertex->x = ds->metadata.hist.mbr.MinX;
    first_Vertex->y = ds->metadata.hist.mbr.MinY;*/
  eulerskew_calculate_bucket_with_mbr(&ds->metadata.hist, first_bucket);
  int buckets = 1;
  listaEulerskew->bucketsList = g_list_append(listaEulerskew->bucketsList, first_bucket);
  /*	listaEulerskew->EdgesList = g_list_append(listaEulerskew->EdgesList, first_Edge);
    listaEulerskew->VertexesList = g_list_append(listaEulerskew->VertexesList, first_Vertex);*/
  double global_skew = first_bucket->skew;
  while (buckets < buckets_num)
  {
    GList *chosen = NULL;
    double skew_reduction = 0;
    GList *item;
    g_list_foreach(item, listaEulerskew->bucketsList)
    {
      eulerskew_face *bucket = (eulerskew_face *)item->data;
      /*		eulerskew_edge *bucketEdge = (eulerskew_edge *)item->data;
          eulerskew_vertex *bucketVertex = (eulerskew_vertex *)item->data;
                bucketEdge->mbr.MinX = bucket->mbr.MinX;
                bucketEdge->mbr.MinY = bucket->mbr.MinY;
                bucketEdge->mbr.MaxX = bucket->mbr.MaxX;
                bucketEdge->mbr.MaxY = bucket->mbr.MaxY;
                bucketVertex->x = bucket->mbr.MinX;
                bucketVertex->y = bucket->mbr.MinY;
      */
      if (isnan(bucket->skew_reduction))
      {
        eulerskew_calculate_skew_reduction(dh, bucket);
      }
      if (skew_reduction < bucket->skew_reduction)
      {
        chosen = item;
        skew_reduction = bucket->skew_reduction;
      }
    }
    // split the chosen bucket
    eulerskew_face *chosen_bucket = (eulerskew_face *)chosen->data;
    eulerskew_face newb1, newb2;
    newb1.mbr = chosen_bucket->mbr;
    newb2.mbr = chosen_bucket->mbr;
    if (chosen_bucket->split_axis == 'x')
      newb1.mbr.MaxX = newb2.mbr.MinX = dh->xtics[chosen_bucket->split_point];
    else
      newb1.mbr.MaxY = newb2.mbr.MinY = dh->ytics[chosen_bucket->split_point];
    eulerskew_calculate_bucket_with_mbr(dh, &newb1);
    eulerskew_calculate_bucket_with_mbr(dh, &newb2);
    global_skew -= chosen_bucket->skew;
    global_skew += newb1.skew + newb2.skew;
    // substitute the old bucket data with b2
    *chosen_bucket = newb1;
    // add the new bucket b2
    eulerskew_face *newb = g_new(eulerskew_face, 1);
    *newb = newb2;
    listaEulerskew->bucketsList = g_list_append(listaEulerskew->bucketsList, newb);
    /*		listaEulerskew->EdgesList = g_list_append(listaEulerskew->EdgesList, newb);
        listaEulerskew->VertexesList = g_list_append(listaEulerskew->VertexesList, newb);
    */
    buckets++;
  }
  // O HISTOGRAMA MINSKEW É APENAS UMA LISTA DE BUCKETS EM QUE AS EXTREMIDADES DO BUCKETS
  // SÃO DEFINIDAS POR DOIS PONTOS minX, minY e maxX, maxY.
  // NO HISTOGRAMA MINSKEW O TAMANHO DAS ARESTAS NÃO É ARMAZENADO.
  // O RESTO DO CÓDIGO DESSA FUNÇÃO ARMAZENA AS DIMENSÕES DAS ARESTAS DE ACORDO COM A LOGICA
  // DO HISTOGRAMA DE EULER.
  GList *item;
  g_list_foreach(item, listaEulerskew->bucketsList)
  {
    eulerskew_face *bucket = (eulerskew_face *)item->data;
    bool edgeExists = false;
    bool vertexExists = false;
    GList *edgeGlist;
    Envelope EdgeMbr;
    // horizontal edge 1
    edgeExists = false;
    EdgeMbr.MinX = bucket->mbr.MinX;
    EdgeMbr.MinY = bucket->mbr.MinY;
    EdgeMbr.MaxX = bucket->mbr.MaxX;
    EdgeMbr.MaxY = bucket->mbr.MinY + 1e-10;
    g_list_foreach(edgeGlist, listaEulerskew->EdgesList)
    {
      eulerskew_edge *edge = (eulerskew_edge *)edgeGlist->data;
      if (ENVELOPE_CONTAINS(EdgeMbr, edge->mbr) || ENVELOPE_CONTAINS(edge->mbr, EdgeMbr))
      {
        edgeExists = true;
        break;
      }
    }
    if (!edgeExists)
    {
      eulerskew_edge *edgeh1 = g_new0(eulerskew_edge, 1);
      edgeh1->mbr = EdgeMbr;
      listaEulerskew->EdgesList = g_list_append(listaEulerskew->EdgesList, edgeh1);
    }
    edgeExists = false;
    // horizontal edge 2
    EdgeMbr.MinX = bucket->mbr.MinX;
    EdgeMbr.MinY = bucket->mbr.MaxY;
    EdgeMbr.MaxX = bucket->mbr.MaxX;
    EdgeMbr.MaxY = bucket->mbr.MaxY + 1e-10;
    g_list_foreach(edgeGlist, listaEulerskew->EdgesList)
    {
      eulerskew_edge *edge = (eulerskew_edge *)edgeGlist->data;
      if (ENVELOPE_CONTAINS(EdgeMbr, edge->mbr) || ENVELOPE_CONTAINS(edge->mbr, EdgeMbr))
      {
        edgeExists = true;
        break;
      }
    }
    if (!edgeExists)
    {
      eulerskew_edge *edgeh2 = g_new0(eulerskew_edge, 1);
      edgeh2->mbr = EdgeMbr;
      listaEulerskew->EdgesList = g_list_append(listaEulerskew->EdgesList, edgeh2);
    }
    edgeExists = false;
    // vertical edge 1
    EdgeMbr.MinX = bucket->mbr.MinX;
    EdgeMbr.MinY = bucket->mbr.MinY;
    EdgeMbr.MaxX = bucket->mbr.MinX + 1e-10;
    EdgeMbr.MaxY = bucket->mbr.MaxY;
    g_list_foreach(edgeGlist, listaEulerskew->EdgesList)
    {
      eulerskew_edge *edge = (eulerskew_edge *)edgeGlist->data;
      if (ENVELOPE_CONTAINS(EdgeMbr, edge->mbr) || ENVELOPE_CONTAINS(edge->mbr, EdgeMbr))
      {
        edgeExists = true;
        break;
      }
    }
    if (!edgeExists)
    {
      eulerskew_edge *edgev1 = g_new0(eulerskew_edge, 1);
      edgev1->mbr = EdgeMbr;
      listaEulerskew->EdgesList = g_list_append(listaEulerskew->EdgesList, edgev1);
    }
    edgeExists = false;
    // vertical edge 2
    EdgeMbr.MinX = bucket->mbr.MaxX;
    EdgeMbr.MinY = bucket->mbr.MinY;
    EdgeMbr.MaxX = bucket->mbr.MaxX + 1e-10;
    EdgeMbr.MaxY = bucket->mbr.MaxY;
    g_list_foreach(edgeGlist, listaEulerskew->EdgesList)
    {
      eulerskew_edge *edge = (eulerskew_edge *)edgeGlist->data;
      if (ENVELOPE_CONTAINS(EdgeMbr, edge->mbr) || ENVELOPE_CONTAINS(edge->mbr, EdgeMbr))
      {
        edgeExists = true;
        break;
      }
    }
    if (!edgeExists)
    {
      eulerskew_edge *edgev2 = g_new0(eulerskew_edge, 1);
      edgev2->mbr = EdgeMbr;
      listaEulerskew->EdgesList = g_list_append(listaEulerskew->EdgesList, edgev2);
    }
    // below left
  

    eulerskew_vertex *vertexbl = g_new0(eulerskew_vertex, 1);
    vertexbl->x = bucket->mbr.MinX;
    vertexbl->y = bucket->mbr.MinY;
    g_list_foreach(edgeGlist, listaEulerskew->VertexesList)
    {
      eulerskew_vertex *vertex = (eulerskew_vertex *)edgeGlist->data;
      if(vertex->x == vertexbl->x && vertex->y == vertexbl->y){
        vertexExists = true;
      }
    }
    if(!vertexExists){
      listaEulerskew->VertexesList = g_list_append(listaEulerskew->VertexesList, vertexbl);
    }
    vertexExists = false;
    // above left
    eulerskew_vertex *vertexal = g_new0(eulerskew_vertex, 1);
    vertexal->x = bucket->mbr.MinX;
    vertexal->y = bucket->mbr.MaxY;
     g_list_foreach(edgeGlist, listaEulerskew->VertexesList)
    {
      eulerskew_vertex *vertex = (eulerskew_vertex *)edgeGlist->data;
      if(vertex->x == vertexal->x && vertex->y == vertexal->y){
        vertexExists = true;
      }
    }
     if(!vertexExists){
    listaEulerskew->VertexesList = g_list_append(listaEulerskew->VertexesList, vertexal);
     }
     vertexExists = false;
    // below right
    eulerskew_vertex *vertexbr = g_new0(eulerskew_vertex, 1);
    vertexbr->x = bucket->mbr.MaxX;
    vertexbr->y = bucket->mbr.MinY;
     g_list_foreach(edgeGlist, listaEulerskew->VertexesList)
    {
      eulerskew_vertex *vertex = (eulerskew_vertex *)edgeGlist->data;
      if(vertex->x == vertexbr->x && vertex->y == vertexbr->y){
        vertexExists = true;
      }
    }
    if(!vertexExists){
    listaEulerskew->VertexesList = g_list_append(listaEulerskew->VertexesList, vertexbr);
    }
    vertexExists = false;
    // above right
    eulerskew_vertex *vertexar = g_new0(eulerskew_vertex, 1);
    vertexar->x = bucket->mbr.MaxX;
    vertexar->y = bucket->mbr.MaxY;
     g_list_foreach(edgeGlist, listaEulerskew->VertexesList)
    {
      eulerskew_vertex *vertex = (eulerskew_vertex *)edgeGlist->data;
      if(vertex->x == vertexar->x && vertex->y == vertexar->y){
        vertexExists = true;
      }
    }
    if(!vertexExists){
    listaEulerskew->VertexesList = g_list_append(listaEulerskew->VertexesList, vertexar);
    }
    vertexExists = false;
  }
  return listaEulerskew;
}
/*
As duas funções a baixo não foram adicionadas em eulerskew.h

double minskew_search_hist(GList *hist, Envelope query) {
  double result = 0.0;
  GList *item;
  g_list_foreach(item, hist) {
    eulerskew_face *b = (eulerskew_face*)item->data;
    if (ENVELOPE_INTERSECTS_EULERSKEW(query, b->mbr)) {
      Envelope inters = EnvelopeIntersection(query, b->mbr);
      double int_area = ENVELOPE_AREA(inters);
      double bucket_area = ENVELOPE_AREA(b->mbr);
      double fraction = int_area / bucket_area;
      result += fraction * b->cardin;
    }
  }
  return result;
}

void eulerskew_print_hist(dataset *ds, GList *hist) {
  char filename[100];
  char *prefix = getenv("HISTOPREFIX");
  prefix = prefix != NULL ? prefix : "";
  sprintf(filename, "histogram/%sminskew-%s.geojson", prefix, ds->metadata.name);
  FILE *f = fopen(filename, "wb");
  if (f == NULL) {
    perror("Error printing histogram");
    return;
  }
  fprintf(f, "{'type': 'FeatureCollection', 'features': [\n");
  int i = 0;
  GList *item;
  g_list_foreach(item, hist) {
    eulerskew_face *bucket = (eulerskew_face *)item->data;
    Envelope e = bucket->mbr;
    fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Polygon\", \"coordinates\": [[");
    fprintf(f, "[%f, %f],", e.MinX, e.MinY);
    fprintf(f, "[%f, %f],", e.MaxX, e.MinY);
    fprintf(f, "[%f, %f],", e.MaxX, e.MaxY);
    fprintf(f, "[%f, %f],", e.MinX, e.MaxY);
    fprintf(f, "[%f, %f]",  e.MinX, e.MinY);
    fprintf(f, "]]}, 'properties': {");
    fprintf(f, "\"name\": \"%d\",", i);
    fprintf(f, "\"skew\": %f,", bucket->skew);
    fprintf(f, "\"card\": %f,", bucket->cardin);
    fprintf(f, "}},\n");
    i++;
  }
  fprintf(f, "]}\n");
  fclose(f);
}
*/
void eulerskew_print_hist(dataset *ds, minskewLists *eh)
{
  char filename[100];
  char *prefix = getenv("HISTOPREFIX");
  prefix = prefix != NULL ? prefix : "";
  sprintf(filename, "histogram/%seulerskew-%s.geojson", prefix, ds->metadata.name);
  FILE *f = fopen(filename, "wb");
  if (f == NULL)
  {
    perror("Error printing histogram");
    return;
  }
  fprintf(f, "{'type': 'FeatureCollection', 'features': [\n");
  GList *item;
  // face
  g_list_foreach(item, eh->bucketsList)
  {
    eulerskew_face *ef = (eulerskew_face *)item->data;
    fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Polygon\", \"coordinates\": [[");
    fprintf(f, "[%.15lf, %.15lf],", ef->mbr.MinX, ef->mbr.MinY);
    fprintf(f, "[%.15lf, %.15lf],", ef->mbr.MaxX, ef->mbr.MinY);
    fprintf(f, "[%.15lf, %.15lf],", ef->mbr.MaxX, ef->mbr.MaxY);
    fprintf(f, "[%.15lf, %.15lf],", ef->mbr.MinX, ef->mbr.MaxY);
    fprintf(f, "[%.15lf, %.15lf]", ef->mbr.MinX, ef->mbr.MinY);
    fprintf(f, "]]}, 'properties': {");
    // fprintf(f, "\"name\": \"f:%d.%d\",", x, y);
    fprintf(f, "\"card\": %lf,", ef->cardin);
    fprintf(f, "\"avg_heigth\": %lf,", ef->avg_height);
    fprintf(f, "\"avg_width\": %lf,", ef->avg_width);
    fprintf(f, "\"avg_area\": %lf,", ef->avg_area);
    fprintf(f, "\"face_area\": %lf,", (ef->mbr.MinX - ef->mbr.MaxX) * (ef->mbr.MinY - ef->mbr.MaxY));
    fprintf(f, "\"type\": \"face\",");
    fprintf(f, "}},\n");
  }
  GList *item2;
  // edges
  // horizontal 1
  g_list_foreach(item2, eh->EdgesList)
  {
    eulerskew_edge *ee = (eulerskew_edge *)item2->data;
    fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"LineString\", \"coordinates\": [[%.15lf, %.15lf], [%.15lf, %.15lf]]},",
            ee->mbr.MinX, ee->mbr.MinY, ee->mbr.MaxX, ee->mbr.MaxY);
    fprintf(f, "'properties': {");
    // fprintf(f, "\"name\": \"eh:%d.%d\",", x, y);
    fprintf(f, "\"card\": %lf,", ee->cardin);
    fprintf(f, "\"type\": \"edge\",");
    fprintf(f, "}},\n");
  }
  // vertex
  g_list_foreach(item, eh->VertexesList)
  {
    eulerskew_vertex *v = (eulerskew_vertex *)item->data;
    fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Point\", \"coordinates\": [%.15lf, %.15lf]},", v->x, v->y);
    fprintf(f, "'properties': {");
    // fprintf(f, "\"name\": \"v:%d.%d\",", v->x, v->y);
    fprintf(f, "\"card\": %lf,", v->cardin);
    fprintf(f, "\"type\": \"vertex\",");
    fprintf(f, "}},\n");
  }
  fprintf(f, "]}\n");
  fclose(f);
}
//
void eulerskew_alloc(dataset *ds, eulerskew_histogram *eh, int xqtd, int yqtd, double psizex, double psizey)
{
  assert(xqtd > 0 && yqtd > 0 && "X and Y must be greater than zero.");
  eh->xqtd = xqtd;
  eh->yqtd = yqtd;
  eh->xtics = g_new(double, eh->xqtd + 1);
  eh->ytics = g_new(double, eh->yqtd + 1);
  eh->faces = g_new0(eulerskew_face, eh->xqtd * eh->yqtd);
  eh->edges = g_new0(eulerskew_edge, ((xqtd + 1) * yqtd) + ((yqtd + 1) * xqtd));
  eh->vertexes = g_new0(eulerskew_vertex, (xqtd + 1) * (yqtd + 1));
  eh->xsize = psizex;
  eh->ysize = psizey;
  printf("Generating histogram of size: %d x %d\n", eh->xqtd, eh->yqtd);
  // X tics
  double xini = ds->metadata.hist.mbr.MinX;
  for (int i = 0; i < eh->xqtd; i++)
    eh->xtics[i] = xini + (psizex * i);
  eh->xtics[eh->xqtd] = ds->metadata.hist.mbr.MaxX;
  // Y tics
  double yini = ds->metadata.hist.mbr.MinY;
  for (int i = 0; i < eh->yqtd; i++)
    eh->ytics[i] = yini + (psizey * i);
  eh->ytics[eh->yqtd] = ds->metadata.hist.mbr.MaxY;
  /*
  GList *item;
  g_list_foreach(item, minskewhist) {
    eulerskew_face *bucket = (eulerskew_face *)item->data;
    //percorre os buckets

    // edges and vertexes
  int v = 0;
  for(int i = 0; i <= eh->xqtd; i++) { // x
      for(int j = 0; j <= eh->yqtd; j++) { // y

          //eulerskew face não tem mbr
          // set vetex at (i,j)
          eh->vertexes[v].x = bucket->mbr.MinX;
          eh->vertexes[v].y = bucket->mbr.MinY;
          v++;
          //os vertices são duas variaveis double, x e y.
          //para cada bucket, o vertice sera o ponto inferior da esquerda.

          // horizontal edge at vertex v
          if (i < eh->xqtd) {
              int e = GET_HORZ_EDGE(i, j);
              eh->edges[e].mbr.MinX = bucket->mbr.MinX;
              eh->edges[e].mbr.MaxX = bucket->mbr.MaxX;
              eh->edges[e].mbr.MinY = bucket->mbr.MinY;
              eh->edges[e].mbr.MaxY = bucket->mbr.MaxY;
          }

          // vertical edge at vertex v
          if (j < eh->yqtd) {
              int e = GET_VERT_EDGE(i, j);
              eh->edges[e].mbr.MinY = bucket->mbr.MinY;
              eh->edges[e].mbr.MaxY = bucket->mbr.MaxY;
              eh->edges[e].mbr.MinX = bucket->mbr.MinX;
              eh->edges[e].mbr.MaxX = bucket->mbr.MaxX;
          }
      }
  }
  }
  */
}
void eulerskew_hash_ds_objects(dataset *ds, eulerskew_histogram *eh, enum JoinPredicateCheck pcheck, minskewLists *ml)
{
  dataset_iter di;
  dataset_foreach(di, ds)
  {
    dataset_leaf *l = get_join_pair_leaf(di.item, pcheck);
    Envelope rs = l->mbr; // pega mbr do objeto
    GEOSGeometryH geo = dataset_get_leaf_geo(ds, l);

    // face
    GList *item;
    g_list_foreach(item, ml->bucketsList)
    {
      eulerskew_face *face = (eulerskew_face *)item->data;
      if (ENVELOPE_CONTAINS(face->mbr, rs) && ENVELOPE_CONTAINS(rs, face->mbr))
      {
        GEOSGeometryH clipped = GEOSClipByRect(geo, face->mbr.MinX, face->mbr.MinY, face->mbr.MaxX, face->mbr.MaxY);
        if (clipped == NULL)
          continue;
        Envelope ev2;
        GEOSEnvelopeGetXY(clipped, &ev2.MinX, &ev2.MaxX, &ev2.MinY, &ev2.MaxY);
        GEOSGeom_destroy(clipped);
        // eulerskew_face *face = &ml->faces[x*ml->yqtd +y];
        // eulerskew_face *face = (eulerskew_face *)item->data;
        face->cardin += 1;
        double delta_x = (ev2.MaxX - ev2.MinX);
        double delta_y = (ev2.MaxY - ev2.MinY);
        double area = delta_x * delta_y;
        face->avg_width += (delta_x - face->avg_width) / face->cardin;
        face->avg_height += (delta_y - face->avg_height) / face->cardin;
        face->avg_area += (area - face->avg_area) / face->cardin;
      }
    }
    // GList *item2;
    //  vertex
    g_list_foreach(item, ml->VertexesList)
    {
      eulerskew_vertex *vertex = (eulerskew_vertex *)item->data;
      // int v = x * (ml->yqtd+1) + y;
      if (ENVELOPE_CONTAINSP(rs, vertex->x, vertex->y))
      {
        vertex->cardin += 1;
      }
    }
    // horizontal edge
    // GList *item3;

    g_list_foreach(item, ml->EdgesList)
    {

      eulerskew_edge *ee = (eulerskew_edge *)item->data;
      if ((ee->mbr.MaxX - ee->mbr.MinX) > (ee->mbr.MaxY - ee->mbr.MinY))
      { // entra no if se a aresta for horizontal
        if (ENVELOPE_INTERSECTS(ee->mbr, rs))
        {
          double delta_x = rs.MaxX - rs.MinX;
          ee->cardin += 1;
          double edge_size = (ee->mbr.MaxX - ee->mbr.MinX);
          ee->avg_projection += (delta_x - ee->avg_projection) / ee->cardin;
        }
      }
      else
      {
        if (ENVELOPE_INTERSECTS(ee->mbr, rs))
        {
          double delta_y = rs.MaxY - rs.MinY;
          ee->cardin += 1;
          double edge_size = (ee->mbr.MaxY - ee->mbr.MinY);
          ee->avg_projection += (delta_y - ee->avg_projection) / ee->cardin;
        }
      }
      // int e = GET_HORZ_EDGE(x, y);
      // ml->avg_projection +=
    }
    // vertical edge

    if (l->gid != -1) // free due to the call to dataset_get_leaf_geo
      GEOSGeom_destroy(geo);
  }
}

void eulerskew_generate_hw(dataset *ds, eulerskew_histogram *eh, double x, double y, enum JoinPredicateCheck pcheck, GList *minskewLists)
{
  double rangex = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
  double rangey = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;
  double psizex = x;
  double psizey = y;
  const int MAX = 1000;
  if ((rangex / psizex) > MAX)
    psizex = rangex / MAX;
  if ((rangey / psizey) > MAX)
    psizey = rangey / MAX;
  eulerskew_alloc(ds, eh, ceil(rangex / psizex), ceil(rangey / psizey), psizex, psizey);
  eulerskew_hash_ds_objects(ds, eh, pcheck, minskewLists);
};
void eulerskew_generate_fix(dataset *ds, eulerskew_histogram *eh, int fsizex, int fsizey, enum JoinPredicateCheck pcheck, GList *minskewLists)
{
  double rangex = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
  double rangey = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;
  double psizex = rangex / fsizex;
  double psizey = rangey / fsizey;
  eulerskew_alloc(ds, eh, fsizex, fsizey, psizex, psizey);
  eulerskew_hash_ds_objects(ds, eh, pcheck, minskewLists);
};
eulerskew_histogram *eulerskew_generate_hist_with_euler(dataset *ds, HistogramGenerateSpec spec, enum JoinPredicateCheck pcheck, GList *minskewLists)
{
  eulerskew_histogram *eh = g_new(eulerskew_histogram, 1);
  eh->mbr = ds->metadata.hist.mbr;
  if (spec.sm == HSPLIT_FIX)
    eulerskew_generate_fix(ds, eh, spec.xqtd, spec.yqtd, pcheck, minskewLists);
  else if (spec.sm == HSPLIT_AVG)
    eulerskew_generate_hw(ds, eh, ds->metadata.x_average, ds->metadata.y_average, pcheck, minskewLists);
  else if (spec.sm == HSPLIT_AVG_STD)
    eulerskew_generate_hw(ds, eh,
                          ds->metadata.x_average + dataset_meta_stddev(ds->metadata, x),
                          ds->metadata.y_average + dataset_meta_stddev(ds->metadata, y),
                          pcheck, minskewLists);
  else
  {
    fprintf(stderr, "Histogram Split Method not implemented.\n");
  }
  return eh;
  //	printf("Generated histogram %d x %d, %s.\n", ds->metadata.hist.xqtd,
  //		ds->metadata.hist.yqtd, HistogramHashMethodName[spec.hm]);
}
int eulerskew_search_hist(eulerskew_histogram *eh, Envelope query2,  minskewLists *listaEulerskew)
{
  if (!ENVELOPE_INTERSECTS(query2, eh->mbr))
    return 0;
  double result = 0;
  printf("result inicio: %f \n", result);
  Envelope query = EnvelopeIntersection2(query2, eh->mbr);
  printf("query mbr: %f %f %f %f, eh mbr: %f %f %f %f \n", query.MaxX, query.MaxY, query.MinX, query.MinY, eh->mbr.MaxX,  eh->mbr.MaxY,  eh->mbr.MinX,  eh->mbr.MinY);
  // face
  GList *item;
  g_list_foreach(item, listaEulerskew->bucketsList)
  {
   
    eulerskew_face *bucket = (eulerskew_face *)item->data;
     printf("cardin %f \n", bucket->cardin);
    printf("result meio2: %f \n", result);
    if (ENVELOPE_CONTAINS(query, bucket->mbr) || ENVELOPE_CONTAINS( bucket->mbr, query))
    {
     
      // eulerskew_face *face = &eh->faces[x * eh->yqtd + y];
      Envelope inters = EnvelopeIntersection(query, bucket->mbr);
      double int_area = ENVELOPE_AREA(inters);
      double face_area = ENVELOPE_AREA(bucket->mbr);
      double fraction = int_area / face_area;
      printf("int_area %f, face_area %f, fraction %f, cardin %f \n", int_area, face_area, fraction, bucket->cardin);
      result += fraction * bucket->cardin;
       printf("result face dentro: %f \n", result);
    }

    g_list_foreach(item, listaEulerskew->EdgesList)
    {
      
      eulerskew_edge *edge = (eulerskew_edge *)item->data;
      eulerskew_edge *edgeNext = (eulerskew_edge *)item->next;
      if ((edge->mbr.MaxX - edge->mbr.MinX) > (edge->mbr.MaxY - edge->mbr.MinY))
      {
        if (ENVELOPE_CONTAINS( query, edge->mbr) || ENVELOPE_CONTAINS( edge->mbr, query))
        {
          //printf("result dentro1: %f \n", result);
          if (edge->mbr.MinY != query.MinY && edgeNext->mbr.MinY != query.MaxY)
          {
            Envelope inters = EnvelopeIntersection(query, edge->mbr);
            // verifica se há intersecção, se há, a função retorna o envelope/mbr da interseção
            double int_length = inters.MaxX - inters.MinX;
            double fraction = int_length / (edge->mbr.MaxX - edge->mbr.MinX);
            result -= fraction * edge->cardin;
            //printf("result dentro2: %f \n", result); o codigo passa por aqui
          }
        }
      }
      else
      {
        if (ENVELOPE_CONTAINS(edge->mbr, query))
        {
          if (edge->mbr.MinX != query.MinX && edgeNext->mbr.MinX != query.MaxX)
          {
            Envelope inters = EnvelopeIntersection(query, edge->mbr);
            double int_length = inters.MaxY - inters.MinY;
            double fraction = int_length / (edge->mbr.MaxY - edge->mbr.MinY);
            result -= fraction * edge->cardin;
          }
        }
      }
    }
  printf("result meio: %f \n", result);
  //   g_list_foreach(item, listaEulerskew->VertexesList)
  //   {
  //     eulerskew_vertex *vertex = (eulerskew_vertex *)item->data;
  //     if (ENVELOPE_CONTAINSP2(query, vertex->x, vertex->y))
  //     {
  //       result += vertex->cardin;
  //     }
  //   }
  // }
  printf("result final: %f \n", result);
return round(result);
}
double eulerskew_estimate_intersections_mp_edges_vert(Envelope el, Envelope er, Envelope inters,
                                                      eulerskew_edge *ehr_face, eulerskew_edge *ehs_face)
{
  // the code below follows equations (1) and (2) in Mamoulis, Papadias 2001
  // estimate the quantity of objects in LeftDs in the inters window: eqn (1)
  double uy = el.MaxY - el.MinY;
  double avgl_y = ehr_face->avg_projection;
  double wy = inters.MaxY - inters.MinY;
  double qtdobjl = ehr_face->cardin *
                   MIN(1, (avgl_y + wy) / uy);
  // estimate the quantity of objects in RightDs in the inters window eqn (1)
  uy = er.MaxY - er.MinY;
  double avgr_y = ehs_face->avg_projection;
  double qtdobjr = ehs_face->cardin *
                   MIN(1, (avgr_y + wy) / uy);
  // estimate join result cardinality, eqn (2)
  return qtdobjl * qtdobjr * MIN(1, (avgl_y + avgr_y) / wy);
}
double eulerskew_estimate_intersections_mp_edges_horz(Envelope el, Envelope er, Envelope inters,
                                                      eulerskew_edge *ehr_face, eulerskew_edge *ehs_face)
{
  // the code below follows equations (1) and (2) in Mamoulis, Papadias 2001
  // estimate the quantity of objects in LeftDs in the inters window: eqn (1)
  double ux = el.MaxX - el.MinX;
  double avgl_x = ehr_face->avg_projection;
  double wx = inters.MaxX - inters.MinX;
  double qtdobjl = ehr_face->cardin *
                   MIN(1, (avgl_x + wx) / ux);
  printf("qtdobjl = %f\n", qtdobjl);
  // estimate the quantity of objects in RightDs in the inters window eqn (1)
  ux = er.MaxX - er.MinX;
  double avgr_x = ehs_face->avg_projection;
  double qtdobjr = ehs_face->cardin *
                   MIN(1, (avgr_x + wx) / ux);
  // estimate join result cardinality, eqn (2)
  return qtdobjl * qtdobjr * MIN(1, (avgl_x + avgr_x) / wx);
}
double eulerskew_estimate_intersections_mamoulis_papadias(Envelope el, Envelope er, Envelope inters,
                                                          eulerskew_face *ehr_face, eulerskew_face *ehs_face)
{
  // the code below follows equations (1) and (2) in Mamoulis, Papadias 2001
  // estimate the quantity of objects in LeftDs in the inters window: eqn (1)
  double ux = el.MaxX - el.MinX;
  double uy = el.MaxY - el.MinY;
  double avgl_x = ehr_face->avg_width;
  double avgl_y = ehr_face->avg_height;
  double wx = inters.MaxX - inters.MinX;
  double wy = inters.MaxY - inters.MinY;
  double qtdobjl = ehr_face->cardin *
                   MIN(1, (avgl_x + wx) / ux) *
                   MIN(1, (avgl_y + wy) / uy);
  // estimate the quantity of objects in RightDs in the inters window eqn (1)
  ux = er.MaxX - er.MinX;
  uy = er.MaxY - er.MinY;
  double avgr_x = ehs_face->avg_width;
  double avgr_y = ehs_face->avg_height;
  double qtdobjr = ehs_face->cardin *
                   MIN(1, (avgr_x + wx) / ux) *
                   MIN(1, (avgr_y + wy) / uy);
  // estimate join result cardinality, eqn (2)
  return qtdobjl * qtdobjr * MIN(1, (avgl_x + avgr_x) / wx) * MIN(1, (avgl_y + avgr_y) / wy);
}
int eulerskew_join_cardinality(dataset *dr,
                               dataset *ds,
                               eulerskew_histogram *ehr,
                               eulerskew_histogram *ehs,
                               rtree_root *rtree_r,
                               rtree_root *rtree_s,
                               double *stddev)
{
  double xini = MAX(ehr->xtics[0], ehs->xtics[0]);
  double yini = MAX(ehr->ytics[0], ehs->ytics[0]);
  double xend = MIN(dr->metadata.hist.mbr.MaxX, ds->metadata.hist.mbr.MaxX);
  double yend = MIN(dr->metadata.hist.mbr.MaxY, ds->metadata.hist.mbr.MaxY);
  unsigned int N = ceil((xend - xini) * (yend - yini));
  double mean = 0.0;
  double M2 = 0.0;
  double sum_error = 0.0;
  int xdr_start = 0;
  while (xdr_start < ehr->xqtd && ehr->xtics[xdr_start + 1] < xini)
    xdr_start++;
  int xdr_end = xdr_start + 1;
  while (xdr_end < ehr->xqtd && ehr->xtics[xdr_end] <= xend)
    xdr_end++;
  if (xdr_start == ehr->xqtd)
    return 0;
  // skip non-intersect area on y
  int ydr_start = 0;
  while (ydr_start < ehr->yqtd && ehr->ytics[ydr_start + 1] < yini)
    ydr_start++;
  int ydr_end = ydr_start + 1;
  while (ydr_end < ehr->yqtd && ehr->ytics[ydr_end] <= yend)
    ydr_end++;
  if (ydr_start == ehr->yqtd)
    return 0;
  int xds_atu = 0;
  float result = 0;
  for (int xr = xdr_start; xr < xdr_end; xr++)
  {
    while (xds_atu < ehs->xqtd && ehs->xtics[xds_atu + 1] < ehr->xtics[xr]) // skip when end of s < start of r
      xds_atu++;
    int xds_end = xds_atu + 1;
    while (xds_end < ehs->xqtd && ehs->xtics[xds_end] <= ehr->xtics[xr + 1]) // increment when end of s < start of r
      xds_end++;
    int yds_atu = 0;
    Envelope er;
    Envelope es;
    er.MinX = ehr->xtics[xr];
    er.MaxX = ehr->xtics[xr + 1];
    for (int yr = ydr_start; yr < ydr_end; yr++)
    {
      while (yds_atu < ehs->yqtd && ehs->ytics[yds_atu + 1] < ehr->ytics[yr]) // skip when end of s < start of r
        yds_atu++;
      int yds_end = yds_atu + 1;
      while (yds_end < ehs->yqtd && ehs->ytics[yds_end] <= ehr->ytics[yr + 1]) // increment when end of s < start of r
        yds_end++;
      er.MinY = ehr->ytics[yr];
      er.MaxY = ehr->ytics[yr + 1];
      double erarea = ENVELOPE_AREA(er);
      for (int xs = xds_atu; xs < xds_end; xs++)
      {
        es.MinX = ehs->xtics[xs];
        es.MaxX = ehs->xtics[xs + 1];
        for (int ys = yds_atu; ys < yds_end; ys++)
        {
          double estimated_result = 0;
          es.MinY = ehs->ytics[ys];
          es.MaxY = ehs->ytics[ys + 1];
          char i = ENVELOPE_INTERSECTS_EULERSKEW(er, es);
          assert(i != 0);
          eulerskew_face *ehr_face = &ehr->faces[xr * ehr->yqtd + yr];
          eulerskew_face *ehs_face = &ehs->faces[xs * ehs->yqtd + ys];
          Envelope inters = EnvelopeIntersection2(er, es);
          double int_area = ENVELOPE_AREA(inters);
          printf("ehr_face[%d][%d].card = %f, avg_h = %f, avg_w = %f\n"
                 "ehs_face[%d][%d].card = %f, avg_h = %f, avg_w = %f\n",
                 xr, yr, ehr_face->cardin, ehr_face->avg_height, ehr_face->avg_width, xs, ys, ehs_face->cardin, ehs_face->avg_height, ehs_face->avg_width);
          double intersections = 0;
          double p = 1;
          if (ehr_face->cardin >= 1 || ehs_face->cardin >= 1)
          {
            intersections = eulerskew_estimate_intersections_mamoulis_papadias(er, es, inters, ehr_face, ehs_face);
            if (intersections < 1.0)
              intersections = 0;
          }
          result += intersections;
          estimated_result += intersections;
          // vertice
          int vr = xr * (ehr->yqtd + 1) + yr;
          int vs = xs * (ehs->yqtd + 1) + ys;
          if (ehs->vertexes[vs].x == ehr->vertexes[vr].x && ehs->vertexes[vs].y == ehr->vertexes[vr].y)
          {
            result += ehr->vertexes[vr].cardin * ehs->vertexes[vs].cardin;
            estimated_result += ehr->vertexes[vr].cardin * ehs->vertexes[vs].cardin;
          }
          // aresta horizontal
          int ar = GET_HORZ_EDGE_EHR(xr, yr);
          int as = GET_HORZ_EDGE_EHS(xs, ys);
          if (ENVELOPE_INTERSECTS_EULERSKEW(ehr->edges[ar].mbr, ehs->edges[as].mbr))
          {
            printf("ar = %d, as = %d\n", ar, as);
            Envelope inters = EnvelopeIntersection2(ehr->edges[ar].mbr, ehs->edges[as].mbr);
            double int_length = inters.MaxX - inters.MinX;
            double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxX - ehr->edges[ar].mbr.MinX);
            double fraction_as = int_length / (ehs->edges[as].mbr.MaxX - ehs->edges[as].mbr.MinX);
            double p = MIN(1, ehr->edges[ar].avg_projection + ehs->edges[as].avg_projection);
            double cardin_ar = ehr->edges[ar].cardin;
            double cardin_as = ehr->edges[as].cardin;
            result -= cardin_ar * cardin_as * p;
            estimated_result -= cardin_ar * cardin_as * p;
            // result -= eulerskew_estimate_intersections_mp_edges_horz(ehr->edges[ar].mbr, ehs->edges[as].mbr, inters, &ehr->edges[ar], &ehs->edges[as]);
          }
          ar = GET_VERT_EDGE_EHR(xr, yr);
          as = GET_VERT_EDGE_EHS(xs, ys);
          if (ENVELOPE_INTERSECTS_EULERSKEW(ehr->edges[ar].mbr, ehs->edges[as].mbr))
          {
            Envelope inters = EnvelopeIntersection2(ehr->edges[ar].mbr, ehs->edges[as].mbr);
            double int_length = inters.MaxY - inters.MinY;
            double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxY - ehr->edges[ar].mbr.MinY);
            double fraction_as = int_length / (ehs->edges[as].mbr.MaxY - ehs->edges[as].mbr.MinY);
            double p = MIN(1, ehr->edges[ar].avg_projection + ehs->edges[as].avg_projection);
            printf(" edge p = %f\n", ehr->edges[ar].avg_projection + ehs->edges[as].avg_projection);
            double cardin_ar = ehr->edges[ar].cardin;
            double cardin_as = ehr->edges[as].cardin;
            result -= cardin_ar * cardin_as * p;
            estimated_result -= cardin_ar * cardin_as * p;
            // result -= eulerskew_estimate_intersections_mp_edges_vert(ehr->edges[ar].mbr, ehs->edges[as].mbr, inters, &ehr->edges[ar], &ehs->edges[as]);
          }
          // double real_cardin = real_cardin_eulerskew_histogram_cell(rtree_r, rtree_s, inters);
          double real_cardin = 1;
          int error = abs(estimated_result - real_cardin);
          double delta = error - mean;
          assert(!isnan(delta));
          mean += delta / (double)N;
          assert(!isnan(mean));
          M2 += delta * (error - mean);
          sum_error += error;
        }
      }
    }
  }
  *stddev = sqrt(M2 / (double)N);
  assert(!isnan(*stddev));
  return round(result);
}
int eulerskew_spatial_join(eulerskew_histogram *ehr, eulerskew_histogram *ehs)
{
  if (!ENVELOPE_INTERSECTS_EULERSKEW(ehr->mbr, ehs->mbr))
    return 0;
  double result = 0;
  double xini = MAX(ehr->xtics[0], ehs->xtics[0]);
  double yini = MAX(ehr->ytics[0], ehs->ytics[0]);
  double xfim = MIN(ehr->mbr.MaxX, ehs->mbr.MaxX);
  double yfim = MIN(ehr->mbr.MaxY, ehs->mbr.MaxY);
  for (int xr = xini; xr <= xfim; xr++)
  {
    Envelope er, es;
    er.MinX = ehr->xtics[xr];
    er.MaxX = ehr->xtics[xr + 1];
    for (int yr = yini; yr <= yfim; yr++)
    {
      er.MinY = ehr->xtics[yr];
      er.MaxY = ehr->xtics[yr + 1];
      for (int xs = xini; xs <= xfim; xs++)
      {
        es.MinX = ehs->xtics[xs];
        es.MaxX = ehs->xtics[xs + 1];
        for (int ys = yini; ys <= yfim; ys++)
        {
          es.MinY = ehs->xtics[ys];
          es.MaxY = ehs->xtics[ys + 1];
          if (ENVELOPE_INTERSECTS_EULERSKEW(er, es))
          {
            // face
            eulerskew_face *ehr_face = &ehr->faces[xr * ehs->yqtd + yr];
            eulerskew_face *ehs_face = &ehr->faces[xs * ehs->yqtd + ys];
            Envelope inters = EnvelopeIntersection(er, es);
            double int_area = ENVELOPE_AREA(inters);
            double face_area = ENVELOPE_AREA(es);
            double fraction = int_area / face_area;
            // result += fraction * face->cardin;
            if (ehr_face->avg_height + ehs_face->avg_height >= 1 && ehr_face->avg_width + ehs_face->avg_width >= 1)
              result += fraction * ehs_face->cardin;
            else
            {
              double p = ehr_face->avg_area + ehs_face->avg_area + ehr_face->avg_height * ehs_face->avg_width + ehs_face->avg_height * ehr_face->avg_width;
              result += fraction * ehs_face->cardin * p;
            }
            // vertice
            int vr = xr * (ehr->yqtd + 1) + yr;
            int vs = xs * (ehs->yqtd + 1) + ys;
            if (ENVELOPE_CONTAINSP(es, ehr->vertexes[vr].x, ehr->vertexes[vr].y))
            {
              result += ehr->vertexes[vr].cardin * ehs->vertexes[vs].cardin;
            }
            // aresta horizontal
            int ar = GET_HORZ_EDGE_EHR(xr, yr);
            int as = GET_HORZ_EDGE_EHS(xs, ys);
            if (ENVELOPE_INTERSECTS_EULERSKEW(ehs->edges[as].mbr, er))
            {
              if (ehs->edges[as].mbr.MinY != es.MinY && ehs->edges[as + 1].mbr.MinY != es.MaxY)
              {
                Envelope inters = EnvelopeIntersection(es, ehs->edges[as].mbr);
                double int_length = inters.MaxX - inters.MinX;
                double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxX - ehr->edges[ar].mbr.MinX);
                double fraction_as = int_length / (ehs->edges[as].mbr.MaxX - ehs->edges[as].mbr.MinX);
                double p = MIN(1, fraction_ar * ehs->edges[ar].cardin + fraction_as * ehs->edges[as].cardin);
                result -= ehr->edges[ar].cardin * ehs->edges[as].cardin * p;
              }
            }
            ar = GET_VERT_EDGE_EHR(xr, yr);
            as = GET_VERT_EDGE_EHS(xs, ys);
            if (ENVELOPE_INTERSECTS_EULERSKEW(ehs->edges[as].mbr, er))
            {
              if (ehs->edges[as].mbr.MinX != es.MinX && ehs->edges[as + 1].mbr.MinX != es.MaxX)
              {
                Envelope inters = EnvelopeIntersection(es, ehs->edges[as].mbr);
                double int_length = inters.MaxY - inters.MinY;
                double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxX - ehr->edges[ar].mbr.MinX);
                double fraction_as = int_length / (ehs->edges[as].mbr.MaxX - ehs->edges[as].mbr.MinX);
                double p = MIN(1, fraction_ar * ehs->edges[ar].cardin + fraction_as * ehs->edges[as].cardin);
                double fraction = int_length / (ehs->edges[as].mbr.MaxY - ehr->edges[as].mbr.MinY);
                result -= ehr->edges[ar].cardin * ehs->edges[as].cardin * p;
              }
            }
          }
        }
      }
    }
  }
  return round(result);
}
int eulerskew_cardinality_per_face(dataset *dr,
                                   dataset *ds,
                                   eulerskew_histogram *ehr,
                                   eulerskew_histogram *ehs,
                                   rtree_root *rtree_r,
                                   rtree_root *rtree_s)
{
  rtree_window_stat stats1;
  rtree_window_stat stats2;
  double xini = MAX(ehr->xtics[0], ehs->xtics[0]);
  double yini = MAX(ehr->ytics[0], ehs->ytics[0]);
  double xend = MIN(dr->metadata.hist.mbr.MaxX, ds->metadata.hist.mbr.MaxX);
  double yend = MIN(dr->metadata.hist.mbr.MaxY, ds->metadata.hist.mbr.MaxY);
  int xdr_start = 0;
  while (xdr_start < ehr->xqtd && ehr->xtics[xdr_start + 1] < xini)
    xdr_start++;
  int xdr_end = xdr_start + 1;
  while (xdr_end < ehr->xqtd && ehr->xtics[xdr_end] <= xend)
    xdr_end++;
  if (xdr_start == ehr->xqtd)
    return 0;
  // skip non-intersect area on y
  int ydr_start = 0;
  while (ydr_start < ehr->yqtd && ehr->ytics[ydr_start + 1] < yini)
    ydr_start++;
  int ydr_end = ydr_start + 1;
  while (ydr_end < ehr->yqtd && ehr->ytics[ydr_end] <= yend)
    ydr_end++;
  if (ydr_start == ehr->yqtd)
    return 0;
  int xds_atu = 0;
  float result = 0;
  for (int xr = xdr_start; xr < xdr_end; xr++)
  {
    while (xds_atu < ehs->xqtd && ehs->xtics[xds_atu + 1] < ehr->xtics[xr]) // skip when end of s < start of r
      xds_atu++;
    int xds_end = xds_atu + 1;
    while (xds_end < ehs->xqtd && ehs->xtics[xds_end] <= ehr->xtics[xr + 1]) // increment when end of s < start of r
      xds_end++;
    int yds_atu = 0;
    Envelope er;
    Envelope es;
    er.MinX = ehr->xtics[xr];
    er.MaxX = ehr->xtics[xr + 1];
    for (int yr = ydr_start; yr < ydr_end; yr++)
    {
      while (yds_atu < ehs->yqtd && ehs->ytics[yds_atu + 1] < ehr->ytics[yr]) // skip when end of s < start of r
        yds_atu++;
      int yds_end = yds_atu + 1;
      while (yds_end < ehs->yqtd && ehs->ytics[yds_end] <= ehr->ytics[yr + 1]) // increment when end of s < start of r
        yds_end++;
      er.MinY = ehr->ytics[yr];
      er.MaxY = ehr->ytics[yr + 1];
      double erarea = ENVELOPE_AREA(er);
      for (int xs = xds_atu; xs < xds_end; xs++)
      {
        es.MinX = ehs->xtics[xs];
        es.MaxX = ehs->xtics[xs + 1];
        for (int ys = yds_atu; ys < yds_end; ys++)
        {
          es.MinY = ehs->ytics[ys];
          es.MaxY = ehs->ytics[ys + 1];
          char i = ENVELOPE_INTERSECTS_EULERSKEW(er, es);
          assert(i != 0);
          eulerskew_face *ehr_face = &ehr->faces[xr * ehr->yqtd + yr];
          eulerskew_face *ehs_face = &ehs->faces[xs * ehs->yqtd + ys];
          Envelope inters = EnvelopeIntersection2(er, es);
          double int_area = ENVELOPE_AREA(inters);
          double estimated_cardin = 0;
          double real_cardin = 0;
          if (ehr_face->cardin >= 1 || ehs_face->cardin >= 1)
          {
            estimated_cardin = eulerskew_estimate_intersections_mamoulis_papadias(er, es, inters, ehr_face, ehs_face);
            real_cardin = real_cardin_eulerskew_histogram_cell(rtree_r, rtree_s, inters);
            char wkt_inters[512];
            sprintf(wkt_inters, "POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
                    inters.MinX, inters.MinY,
                    inters.MaxX, inters.MinY,
                    inters.MaxX, inters.MaxY,
                    inters.MinX, inters.MaxY,
                    inters.MinX, inters.MinY);
            GEOSGeometryH geo_inters = GEOSGeomFromWKT(wkt_inters);
            memset(&stats1, 0, sizeof(rtree_window_stat));
            memset(&stats2, 0, sizeof(rtree_window_stat));
            GList *results_ehr = rtree_window_search(rtree_r, geo_inters, &stats1);
            GList *results_ehs = rtree_window_search(rtree_s, geo_inters, &stats2);
            GList *a;
            g_list_foreach(a, results_ehr)
            {
              rtree_leaf *la = (rtree_leaf *)a->data;
              GList *b;
              g_list_foreach(b, results_ehs)
              {
                rtree_leaf *lb = (rtree_leaf *)b->data;
                if (GEOSIntersects(la->geo, lb->geo))
                {
                  /*Envelope eva; GEOSGetEnvelope(la->geo, &eva);
                  Envelope evb; GEOSGetEnvelope(lb->geo, &evb);
                  Envelope inters = EnvelopeIntersection2(eva, evb);
                  double area = ENVELOPE_AREA(inters);*/
                  real_cardin += 1;
                }
              }
            }
            if (estimated_cardin < 1.0)
              estimated_cardin = 0;
          }
          int vr = xr * (ehr->yqtd + 1) + yr;
          int vs = xs * (ehs->yqtd + 1) + ys;
          if (ehs->vertexes[vs].x == ehr->vertexes[vr].x && ehs->vertexes[vs].y == ehr->vertexes[vr].y)
          {
            result += ehr->vertexes[vr].cardin * ehs->vertexes[vs].cardin;
          }
          // aresta horizontal
          int ar = GET_HORZ_EDGE_EHR(xr, yr);
          int as = GET_HORZ_EDGE_EHS(xs, ys);
          if (ENVELOPE_INTERSECTS_EULERSKEW(ehr->edges[ar].mbr, ehs->edges[as].mbr))
          {
            Envelope inters = EnvelopeIntersection2(ehr->edges[ar].mbr, ehs->edges[as].mbr);
            double int_length = inters.MaxX - inters.MinX;
            double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxX - ehr->edges[ar].mbr.MinX);
            double fraction_as = int_length / (ehs->edges[as].mbr.MaxX - ehs->edges[as].mbr.MinX);
            double p = MIN(1, ehr->edges[ar].avg_projection + ehs->edges[as].avg_projection);
            double cardin_ar = ehr->edges[ar].cardin;
            double cardin_as = ehr->edges[as].cardin;
            result -= cardin_ar * cardin_as * p;
            // result -= eulerskew_estimate_intersections_mp_edges_horz(ehr->edges[ar].mbr, ehs->edges[as].mbr, inters, &ehr->edges[ar], &ehs->edges[as]);
          }
          ar = GET_VERT_EDGE_EHR(xr, yr);
          as = GET_VERT_EDGE_EHS(xs, ys);
          if (ENVELOPE_INTERSECTS_EULERSKEW(ehr->edges[ar].mbr, ehs->edges[as].mbr))
          {
            Envelope inters = EnvelopeIntersection2(ehr->edges[ar].mbr, ehs->edges[as].mbr);
            double int_length = inters.MaxY - inters.MinY;
            double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxY - ehr->edges[ar].mbr.MinY);
            double fraction_as = int_length / (ehs->edges[as].mbr.MaxY - ehs->edges[as].mbr.MinY);
            double p = MIN(1, ehr->edges[ar].avg_projection + ehs->edges[as].avg_projection);
            double cardin_ar = ehr->edges[ar].cardin;
            double cardin_as = ehr->edges[as].cardin;
            result -= cardin_ar * cardin_as * p;
            // result -= eulerskew_estimate_intersections_mp_edges_vert(ehr->edges[ar].mbr, ehs->edges[as].mbr, inters, &ehr->edges[ar], &ehs->edges[as]);
          }
          printf("%d \t %d\n", (int)estimated_cardin, (int)real_cardin);
        }
      }
    }
  }
  return round(result);
}
int real_cardin_eulerskew_histogram_cell(rtree_root *rtree_r, rtree_root *rtree_s, Envelope inters)
{
  rtree_window_stat stats1;
  rtree_window_stat stats2;
  double real_cardin = 0;
  char wkt_inters[512];
  sprintf(wkt_inters, "POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
          inters.MinX, inters.MinY,
          inters.MaxX, inters.MinY,
          inters.MaxX, inters.MaxY,
          inters.MinX, inters.MaxY,
          inters.MinX, inters.MinY);
  GEOSGeometryH geo_inters = GEOSGeomFromWKT(wkt_inters);
  memset(&stats1, 0, sizeof(rtree_window_stat));
  memset(&stats2, 0, sizeof(rtree_window_stat));
  GList *results_ehr = rtree_window_search(rtree_r, geo_inters, &stats1);
  GList *results_ehs = rtree_window_search(rtree_s, geo_inters, &stats2);
  GList *a;
  g_list_foreach(a, results_ehr)
  {
    rtree_leaf *la = (rtree_leaf *)a->data;
    GList *b;
    g_list_foreach(b, results_ehs)
    {
      rtree_leaf *lb = (rtree_leaf *)b->data;
      if (GEOSIntersects(la->geo, lb->geo))
        real_cardin++;
    }
  }
  return round(real_cardin);
}
