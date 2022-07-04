#include "../include/SpatialHistogramAB.hpp"

#include <array>
#include <algorithm>

Cell SpatialHistogramAB::getCell(double x, double y)
{
    Cell cell;
    cell.X = std::floor((x - bottomLeft.X) / cellWidth);
    cell.Y = std::floor((y - bottomLeft.Y) / cellHeight);
    return cell;
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

    bottomLeft.X = metadata.mbr.MinX;
    bottomLeft.Y = metadata.mbr.MinY;

    for(const DatasetEntry& object : ds.geoms()) {
        const Envelope& mbr = object.mbr;
        bool found = false;

        Cell minCell = getCell(mbr.MinX, mbr.MinY);
        Cell maxCell = getCell(mbr.MaxX, mbr.MaxY);

        for (ABBucket& bucket : buckets) {
            if (bucket.contains(minCell, maxCell)) {
                bucket.Cardinality++;
                found = true;
            }
        }

        if (!found) {
            ABBucket newBucket;
            newBucket.MinCell = minCell;
            newBucket.MaxCell = maxCell;
            newBucket.Cardinality = 1;
            newBucket.ID = buckets.size();
            buckets.push_back(newBucket);
        }
    }
}

SpatialHistogramAB::~SpatialHistogramAB()
{
    buckets.clear();
}

Pos SpatialHistogramAB::getCellBottomLeft(Cell cell)
{
    Pos pos;
    pos.X = bottomLeft.X + cell.X * cellWidth;
    pos.Y = bottomLeft.Y + cell.Y * cellHeight;
    return pos;
}

Pos SpatialHistogramAB::getCellTopRight(Cell cell)
{
    Pos pos;
    pos.X = bottomLeft.X + cell.X * cellWidth + cellWidth;
    pos.Y = bottomLeft.Y + cell.Y * cellHeight + cellHeight;
    return pos;
}

Envelope SpatialHistogramAB::getOuterRect(ABBucket bucket)
{
    Pos minPos = getCellBottomLeft(bucket.MinCell);
    Pos maxPos = getCellTopRight(bucket.MaxCell);

    Envelope outer;
    outer.MinX = minPos.X;
    outer.MinY = minPos.Y;
    outer.MaxX = maxPos.X;
    outer.MaxY = maxPos.Y;
    return outer;
}

Envelope SpatialHistogramAB::getInnerRect(ABBucket bucket)
{
    Envelope inner;

    Pos min = getCellBottomLeft(bucket.MinCell);
    Pos max = getCellTopRight(bucket.MaxCell);

    if ((bucket.MaxCell.X - bucket.MinCell.X <= 1) ||
        (bucket.MaxCell.Y - bucket.MinCell.Y <= 1))
    {
        double middleX = min.X + (abs(max.X-min.X)/2.0);
        double middleY = min.Y + (abs(max.Y-min.Y)/2.0);

        inner.MinX = middleX;
        inner.MaxX = middleX;
        inner.MinY = middleY;
        inner.MaxY = middleY;
    }
    else
    {
        inner.MinX = min.X;
        inner.MinY = min.Y;
        inner.MaxX = max.X;
        inner.MaxY = max.Y;
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