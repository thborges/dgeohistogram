/*
 * Envelope structure to represent MBRs
 *
 *  Created on: 2022-02-15
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include <numeric>

struct Envelope { 
	double MinX;
	double MinY;
	double MaxX;
	double MaxY;

	Envelope() {
		MinX = MinY = std::numeric_limits<double>::max();
		MaxX = MaxY = -std::numeric_limits<double>::max();
	}

	Envelope(double minX, double minY, double maxX, double maxY): MinX(minX), MinY(minY), MaxX(maxX), MaxY(maxY) {}

	void merge(const Envelope& e) {
		if (isEmpty())
			*this = e;
		else {
			this->MinX = std::min(this->MinX, e.MinX);
			this->MinY = std::min(this->MinY, e.MinY);
			this->MaxX = std::max(this->MaxX, e.MaxX);
			this->MaxY = std::max(this->MaxY, e.MaxY);
		}
	}

	void merge(double MinX, double MinY, double MaxX, double MaxY) {
		if (isEmpty()) {
			this->MinX = MinX;
			this->MinY = MinY;
			this->MaxX = MaxX;
			this->MaxY = MaxY;
		} else {
			this->MinX = std::min(this->MinX, MinX);
			this->MinY = std::min(this->MinY, MinY);
			this->MaxX = std::max(this->MaxX, MaxX);
			this->MaxY = std::max(this->MaxY, MaxY);
		}
	}

	bool intersects(const Envelope& mbr2) const {
		return !isEmpty() && !mbr2.isEmpty() &&
			MinX <= mbr2.MaxX && MinY <= mbr2.MaxY && mbr2.MinX <= MaxX && mbr2.MinY <= MaxY;
	}

	bool contains(double x, double y) const {
		return !isEmpty() && x >= MinX && x <= MaxX && y >= MinY && y <= MaxY;
	}

	bool contains(const Envelope& mbr2) const {
		return !isEmpty() && !mbr2.isEmpty()
			   && contains(mbr2.MinX, mbr2.MinY)
			   && contains(mbr2.MinX, mbr2.MaxY)
			   && contains(mbr2.MaxX, mbr2.MinY)
			   && contains(mbr2.MaxX, mbr2.MaxY);
	}

	Envelope intersection(const Envelope& r) const {
		if (intersects(r)) {
			return Envelope{
				r.MinX > MinX ? r.MinX : MinX,
				r.MinY > MinY ? r.MinY : MinY,
				r.MaxX < MaxX ? r.MaxX : MaxX,
				r.MaxY < MaxY ? r.MaxY : MaxY};
		}
		else
			return Envelope();
	}

	bool inline isEmpty() const {
		return MaxX < MinX || MaxY < MinY;
	}

	double area() const {
		if (isEmpty())
			return 0.0;
		else
			return (MaxX - MinX) * (MaxY - MinY);
	}

	double width() const {
		if (isEmpty())
			return 0.0;
		else
			return MaxX - MinX;
	}

	double length() const {
		if (isEmpty())
			return 0.0;
		else
			return MaxY - MinY;
	}
};
