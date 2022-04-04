#include "../include/SpatialHistogramAB.hpp"

#include <algorithm>

ABBucket SpatialHistogramAB::getMBRBucket(const Envelope& mbr)
{
    for (ABBucket& bucket : buckets) {
        if (bucket.contains(mbr)) {
            bucket.cardin++;
            return bucket;
        }
    }

    ABBucket newBucket;
    newBucket.cardin = 0;

    int minCellX = std::floor((mbr.MinX - bottomLeftX) / cellWidth);
    int minCellY = std::floor((mbr.MinY - bottomLeftY) / cellHeight);

    int maxCellX = std::floor((mbr.MaxX - bottomLeftX) / cellWidth);
    int maxCellY = std::floor((mbr.MaxY - bottomLeftY) / cellHeight);

    newBucket.outer.MinX = bottomLeftX + minCellX * cellWidth;
    newBucket.outer.MinY = bottomLeftY + minCellY * cellHeight;

    if (minCellX == maxCellX && minCellY == maxCellY)
    {
        newBucket.outer.MaxX = newBucket.outer.MinX + cellWidth;
        newBucket.outer.MaxY = newBucket.outer.MinY + cellHeight;

        const double centerX = newBucket.outer.MinX + (cellWidth/2.0);
        const double centerY = newBucket.outer.MinY + (cellHeight/2.0);

        newBucket.inner.MinX = centerX;
        newBucket.inner.MinY = centerY;
        newBucket.inner.MaxX = centerX;
        newBucket.inner.MaxY = centerY;
    }
    else if ((maxCellX - minCellX) == 1 || (maxCellY - minCellY) == 1)
    {
        newBucket.outer.MaxX = bottomLeftX + maxCellX * cellWidth + cellWidth;
        newBucket.outer.MaxY = bottomLeftY + maxCellY * cellHeight + cellHeight;

        const double centerX = (newBucket.outer.MaxX - newBucket.outer.MinX) / 2.0;
        const double centerY = (newBucket.outer.MaxY - newBucket.outer.MinY) / 2.0;

        newBucket.inner.MinX = centerX;
        newBucket.inner.MinY = centerY;
        newBucket.inner.MaxX = centerX;
        newBucket.inner.MaxY = centerY;
    }
    else
    {
        newBucket.outer.MaxX = bottomLeftX + maxCellX * cellWidth + cellWidth;
        newBucket.outer.MaxY = bottomLeftY + maxCellY * cellHeight + cellHeight;

        newBucket.inner.MinX = newBucket.outer.MinX + cellWidth;
        newBucket.inner.MinY = newBucket.outer.MinY + cellHeight;
        newBucket.inner.MaxX = newBucket.outer.MaxX - cellWidth;
        newBucket.inner.MaxY = newBucket.outer.MaxY - cellHeight;
    }

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
        return a.mbr.area() < b.mbr.area();
    });

    for(const DatasetEntry& object : objects) {
        ABBucket bucket = getMBRBucket(object.mbr);
        if (bucket.cardin == 0) {
            bucket.id = buckets.size();
            bucket.cardin++;
            buckets.push_back(bucket);
        }
    }
}

SpatialHistogramAB::~SpatialHistogramAB()
{
    buckets.clear();
}

double SpatialHistogramAB::intersectionEstimation(const Envelope& wquery)
{
    double result = 0.0;

    for (ABBucket& bucket : buckets) {
        double probability;

        if (!wquery.intersects(bucket.outer)) {
            probability = 0.0;
        } else if (!wquery.intersects(bucket.inner)) {
            Envelope intersection = wquery.intersection(bucket.outer);
            probability = std::min<double>(intersection.width()/cellWidth, intersection.length()/cellHeight);
            if (probability > 1.0)
                probability = 1.0;
            assert(probability >= 0.0 && probability <= 1.0);
        } else {
            probability = 1.0;
        }

        result += probability * bucket.cardin;
    }

    return result;
}

double SpatialHistogramAB::withinEstimation(const Envelope& wquery)
{
    double result = 0.0;

    for (ABBucket& bucket : buckets) {
        double probability;

        if (!wquery.contains(bucket.inner) && !wquery.contains(bucket.outer)) {
            probability = 0.0;
        } else if (!wquery.contains(bucket.outer)) {
            std::array<double, 4> gaps;
            gaps[0] = std::abs(wquery.MaxY - bucket.inner.MaxY) * cellHeight;
            gaps[1] = std::abs(wquery.MaxX - bucket.inner.MaxX) * cellWidth;
            gaps[2] = std::abs(wquery.MinY - bucket.inner.MinY) * cellHeight;
            gaps[3] = std::abs(wquery.MinX - bucket.inner.MinX) * cellWidth;
            probability = *std::min_element(std::begin(gaps), std::end(gaps));
        } else {
            probability = 1.0;
        }

        result += probability * bucket.cardin;
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
        if (!first)
            output << ",\n";
        first = false;
        output << "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Polygon\", \"coordinates\": [";
        output << "[[" << bucket.outer.MinX << "," << bucket.outer.MinY << "],";
        output << "["  << bucket.outer.MaxX << "," << bucket.outer.MinY << "],";
        output << "["  << bucket.outer.MaxX << "," << bucket.outer.MaxY << "],";
        output << "["  << bucket.outer.MinX << "," << bucket.outer.MaxY << "],";
        output << "["  << bucket.outer.MinX << "," << bucket.outer.MinY << "]],";
        output << "[[" << bucket.inner.MinX << "," << bucket.inner.MinY << "],";
        output << "["  << bucket.inner.MaxX << "," << bucket.inner.MinY << "],";
        output << "["  << bucket.inner.MaxX << "," << bucket.inner.MaxY << "],";
        output << "["  << bucket.inner.MinX << "," << bucket.inner.MaxY << "],";
        output << "["  << bucket.inner.MinX << "," << bucket.inner.MinY << "]]";
        output << "]}, \"properties\": {";
        output << "\"name\": \"" << bucket.id << "\",";
        output << "\"card\": " << bucket.cardin << "}}\n";
	}

	output << "]}\n";
	output.close();
}