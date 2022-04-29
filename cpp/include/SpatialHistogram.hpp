/*
 * Spatial Histogram Interface
 *
 *  Created on: 2022-02-15
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include <string>
#include "Envelope.hpp"

class SpatialHistogram {
public:
	virtual ~SpatialHistogram() {};

	virtual void printGeoJson(const std::string& filename) = 0;
	virtual double estimateWQuery(const Envelope& wquery) = 0;
	virtual double getSize() const = 0;
	virtual const std::string name() = 0;

	virtual Envelope const & mbr() { return hmbr; };

protected:
	Envelope hmbr;
};

