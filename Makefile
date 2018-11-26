
INCLUDES= \
	-DHAVE_ISNAN=1 \
	-I/opt/local/include \
	-I/usr/include/gdal \
	-Iauxs

LIBS= \
	-L/opt/local/lib \
	-L/usr/lib

DEBUG=-Wall -ggdb -O0
#DEBUG=-O2

DEPEND_CPP=$(patsubst %.cpp,%.o,$(wildcard auxs/*.cpp)) $(patsubst %.cpp,%.o,$(wildcard *.cpp)) 
DEPEND_C=$(patsubst %.c,%.o,$(wildcard auxs/*.c)) $(patsubst %.c,%.o,$(wildcard *.c)) 

all: main

main: $(DEPEND_CPP) $(DEPEND_C)
	g++ -std=c++11 $(DEBUG) $(INCLUDES) $(DEPEND_C) $(DEPEND_CPP) -o main $(LIBS) -lgdal -lgeos -lgeos_c

%.o: %.c
	gcc -std=c11 $(DEBUG) $(INCLUDES) -c $< -o $@

%.o: %.cpp
	g++ -std=c++11 $(DEBUG) $(INCLUDES) -c $< -o $@

clean:
	rm -f *.o
	rm -f auxs/*.o
	rm -f main

clean_files:
	rm -f mbrs/*.geojson
	rm -f histogram/*.geojson
	rm -f histogram/*.dat
