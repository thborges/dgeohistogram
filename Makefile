
INCLUDES= \
	-D_GNU_SOURCE \
	-I/opt/local/include \
	-I/usr/include/gdal \
	-Iauxs

LIBS= \
	-L/opt/local/lib \
	-L/usr/lib

all: geosext.o
	gcc -std=c11 -ggdb -O0 -c $(INCLUDES) *.c auxs/*.c 
	g++ -std=c++11 -ggdb -O0 *.o -o main $(LIBS) -lgdal -lgeos -lgeos_c


geosext.o: 
	g++ -ggdb -O0 -std=c++11 -c $(INCLUDES) auxs/*.cpp

clean:
	rm -f *.o
	rm -f main

clean_files:
	rm -f mbrs/*.geojson
	rm -f histogram/*.geojson
	rm -f histogram/*.dat