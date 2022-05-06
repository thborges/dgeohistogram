#include "../include/SpatialHistogramAB.hpp"

#include <algorithm>

void SpatialHistogramAB::getCell(double x, double y, int* cellX, int* cellY)
{
    *cellX = std::floor((x - bottomLeftX) / cellWidth);
    *cellY = std::floor((y - bottomLeftY) / cellHeight);
}

ABBucket SpatialHistogramAB::getMBRBucket(const Envelope& mbr)
{
    int minCellX, minCellY;
    getCell(mbr.MinX, mbr.MinY, &minCellX, &minCellY);
    int maxCellX, maxCellY;
    getCell(mbr.MaxX, mbr.MaxY, &maxCellX, &maxCellY);

    for (ABBucket& bucket : buckets) {
        if (bucket.contains(minCellX, minCellY, maxCellX, maxCellY)) {
            bucket.Cardinality++;
            return bucket;
        }
    }

    ABBucket newBucket;
    newBucket.MinCellX = minCellX;
    newBucket.MinCellY = minCellY;
    newBucket.MaxCellX = maxCellX;
    newBucket.MaxCellY = maxCellY;
    newBucket.Cardinality = 0;

    return newBucket;
}

SpatialHistogramAB::SpatialHistogramAB(Dataset& ds, int columns, int rows)
    : columns(columns),
      rows(rows)
{
    const DatasetMetadata& metadata = ds.metadata();

    double rangex = metadata.mbr.MaxX - metadata.mbr.MinX;
    double rangey = metadata.mbr.MaxY - metadata.mbr.MinY;

    cellWidth = rangex / columns;
    cellHeight = rangey / rows;

    bottomLeftX = metadata.mbr.MinX;
    bottomLeftY = metadata.mbr.MinY;

    auto objects = ds.geoms();
    objects.sort([](const DatasetEntry& a, const DatasetEntry& b)
    {
        return a.mbr.area() > b.mbr.area();
    });

    for(const DatasetEntry& object : objects) {
        ABBucket bucket = getMBRBucket(object.mbr);
        if (bucket.Cardinality == 0) {
            bucket.ID = buckets.size();
            bucket.Cardinality++;
            buckets.push_back(bucket);
        }
    }
}

SpatialHistogramAB::~SpatialHistogramAB()
{
    buckets.clear();
}

void SpatialHistogramAB::getCellBottomLeft(int cellX, int cellY, double* x, double* y)
{
    *x = bottomLeftX + cellX * cellWidth;
    *y = bottomLeftY + cellY * cellHeight;
}

void SpatialHistogramAB::getCellTopRight(int cellX, int cellY, double* x, double* y)
{
    *x = bottomLeftX + cellX * cellWidth + cellWidth;
    *y = bottomLeftY + cellY * cellHeight + cellHeight;
}

Envelope SpatialHistogramAB::getOuterRect(ABBucket bucket)
{
    Envelope outer;
    getCellBottomLeft(bucket.MinCellX, bucket.MinCellY, &outer.MinX, &outer.MinY);
    getCellTopRight(bucket.MaxCellX, bucket.MaxCellY, &outer.MaxX, &outer.MaxY);
    return outer;
}

Envelope SpatialHistogramAB::getInnerRect(ABBucket bucket)
{
    Envelope inner;
    if ((bucket.MaxCellX - bucket.MinCellX <= 1) ||
        (bucket.MaxCellY - bucket.MinCellY <= 1))
    {
        double minX, minY;
        double maxX, maxY;
        getCellBottomLeft(bucket.MinCellX, bucket.MinCellY, &minX, &minY);
        getCellTopRight(bucket.MaxCellX, bucket.MaxCellY, &maxX, &maxY);

        double middleX = minX + (abs(maxX-minX)/2.0);
        double middleY = minY + (abs(maxY-minY)/2.0);
        
        inner.MinX = middleX;
        inner.MaxX = middleX;
        inner.MinY = middleY;
        inner.MaxY = middleY;
    }
    else
    {
        getCellTopRight(bucket.MinCellX, bucket.MinCellY, &inner.MinX, &inner.MinY);
        getCellBottomLeft(bucket.MaxCellX, bucket.MaxCellY, &inner.MaxX, &inner.MaxY);
    }
    return inner;
}

double SpatialHistogramAB::intersectionEstimation(const Envelope& wquery)
{
    double result = 0.0;

    for (ABBucket& bucket : buckets) {
        double probability;

        Envelope inner = getInnerRect(bucket);
        Envelope outer = getOuterRect(bucket);

        if (!wquery.intersects(outer)) {
            probability = 0.0;
        } else if (!wquery.intersects(inner)) {
            Envelope intersection = wquery.intersection(outer);
            probability = std::min<double>(intersection.width()/cellWidth, intersection.length()/cellHeight);
            if (probability > 1.0)
                probability = 1.0;
            assert(probability >= 0.0 && probability <= 1.0);
        } else {
            probability = 1.0;
        }

        result += probability * bucket.Cardinality;
    }

    return result;
}

double SpatialHistogramAB::withinEstimation(const Envelope& wquery)
{
    double result = 0.0;

    for (ABBucket& bucket : buckets) {
        double probability;

        Envelope inner = getInnerRect(bucket);
        Envelope outer = getOuterRect(bucket);

        if (!wquery.contains(inner) && !wquery.contains(outer)) {
            probability = 0.0;
        } else if (!wquery.contains(outer)) {
            std::array<double, 4> gaps;
            gaps[0] = std::abs(wquery.MaxY - inner.MaxY) * cellHeight;
            gaps[1] = std::abs(wquery.MaxX - inner.MaxX) * cellWidth;
            gaps[2] = std::abs(wquery.MinY - inner.MinY) * cellHeight;
            gaps[3] = std::abs(wquery.MinX - inner.MinX) * cellWidth;
            probability = *std::min_element(std::begin(gaps), std::end(gaps));
        } else {
            probability = 1.0;
        }

        result += probability * bucket.Cardinality;
    }

    return result;
}

double SpatialHistogramAB::estimateWQuery(const Envelope& wquery)
{
    return intersectionEstimation(wquery);
}

void SpatialHistogramAB::printGeoJson(const std::string& filename)
{
	std::ofstream output;
	output.open(filename);
	output << std::fixed;
	output.precision(15);
	output << "{\"type\": \"FeatureCollection\", \"features\": [\n";

	bool first = true;
	for(ABBucket& bucket : buckets) {
        Envelope inner = getInnerRect(bucket);
        Envelope outer = getOuterRect(bucket);

        if (!first)
            output << ",\n";
        first = false;
        output << "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Polygon\", \"coordinates\": [";
        output << "[[" << outer.MinX << "," << outer.MinY << "],";
        output << "["  << outer.MaxX << "," << outer.MinY << "],";
        output << "["  << outer.MaxX << "," << outer.MaxY << "],";
        output << "["  << outer.MinX << "," << outer.MaxY << "],";
        output << "["  << outer.MinX << "," << outer.MinY << "]],";
        output << "[[" << inner.MinX << "," << inner.MinY << "],";
        output << "["  << inner.MaxX << "," << inner.MinY << "],";
        output << "["  << inner.MaxX << "," << inner.MaxY << "],";
        output << "["  << inner.MinX << "," << inner.MaxY << "],";
        output << "["  << inner.MinX << "," << inner.MinY << "]]";
        output << "]}, \"properties\": {";
        output << "\"name\": \"" << bucket.ID << "\",";
        output << "\"card\": " << bucket.Cardinality << "}}\n";
	}

	output << "]}\n";
	output.close();
}