# README #

This README would normally document whatever steps are necessary to get your application up and running.

### How do I get set up? ###

apt-get install libgdal-dev libgeos++-dev build-essential autotools-dev autoconf
make # to build the sources

You will need at least libgeos++-dev version 3.5. If your linux distribution does not provide it, compile from sources.
Download http://download.osgeo.org/geos/geos-3.5.0.tar.bz2

tar xf geos-3.5.0.tar.bz2

./configure --prefix=/usr/

make

make install (as root)

To configure the dgeohistogram, use, in the main folder:

./configure

make
