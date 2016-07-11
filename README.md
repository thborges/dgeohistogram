# README #

This README would normally document whatever steps are necessary to get your application up and running.

## How do I get set up? ##

### Linux ###

apt-get install libgdal-dev libgeos++-dev build-essential autotools-dev autoconf
make # to build the sources

You will need at least libgeos++-dev version 3.5. If your linux distribution does not provide it, compile from sources.
Download http://download.osgeo.org/geos/geos-3.5.0.tar.bz2

tar xf geos-3.5.0.tar.bz2

./configure --prefix=/usr/

make

make install (as root)

### Windows ###

* Building on Windows x64 ( Tested using Windows 10 x64 )
	
	* We need to use cross compilling because GDAL (2.1.0) doesnt not compile in Windows directly
	* WARNING: Only use CFLAGS and CXXFLAGS with -pipe option if you have enough memory on your environment
	
	* You can build it at your own, or grab the compiled binaries from Dropbox Directory
	
		* Binaries ( Mingw 5.3.x )
			* https://www.dropbox.com/sh/j81byas8y365mkw/AADHh3wGcUyuYQVX78uosf20a?dl=0
	
		* Building (Easy Way - Cross Compilling)
			* OS: Ubuntu Latest >= 16.04 x64 ( Tested using Ubuntu 16.04 Xenial Mini x64 )
		
			* Dependences:
				* sudo apt-get install build-essential gcc-mingw-w64 g++-mingw-w64 mingw-w64 make wget
			
			* Cross Compilling GEOS >= 3.5.0
				* Download source files from: https://trac.osgeo.org/geos/
			
				* Step by Step for version 3.5.0
					* wget http://goo.gl/oWvZW3
					* tar -jxvf oWvZW3
					* cd geos-3.5.0
					* ./configure --host=x86_64-w64-mingw32 CFLAGS="-pipe" CXXFLAGS="-pipe"
					* make -j 8
					* cd ..
					* tar -zcvf GEOS-3.5.0_x64.tar.gz geos-3.5.0
			
			* Cross Compilling GDAL >= 2.1.0
				* Download source files from: https://trac.osgeo.org/gdal/wiki/DownloadSource
			
				* Step by Step for version 2.1.0
					* wget http://goo.gl/ud5zZF
					* tar -zxvf ud5zZF
					* cd gdal-2.1.0
					* ./configure --host=x86_64-w64-mingw32 CFLAGS="-pipe" CXXFLAGS="-pipe"
					* make -j 8
					* cd ..
					* tar -zcvf GDAL-2.1.0_x64.tar.gz gdal-2.1.0
				
			* Now copy GEOS-3.5.0_x64.tar.gz and GDAL-2.1.0_x64.tar.gz to Windows Machine
		
	* Preparing the Workspace
		* Download and Install Mingw64 >= 5.3.0
			* https://sourceforge.net/projects/mingw-w64/
		* Past the binaries (GEOS-3.5.0_x64.tar.gz and GDAL-2.1.0_x64.tar.gz) at folder name "win64" inside dgeohistogram root
		* Extract them there ( Ignore any Symbolic Link Error )
		* Copy all those dlls to bin folder:
			* geos-3.5.0/src/.libs/libgeos-3-5-0.dll
			* geos-3.5.0/capi/.libs/libgeos_c-1.dll
			* gdal-2.1.0/.libs/libgdal-20.dll
			
	* Compilling dgeohistogram
		* This Workspace should work on Linux/Mac with GCC and G++, just ignore the related Windows steps
		* Open CodeLite IDE >= 9.1.8
			* First time the CodeLite IDE bootup it will check the compilers installed in your computer, set Mingw64 to default if you are using Windows OS.
			* Goto Workspace -> Open Workspace
			* Navigate to win64 folder and open the Workspace "dgeohistogram.workspace"
			* Now goto Build -> Build Project
			* Done, now you should have the exe in the bin folder
			