
INCLUDES= \
	-D_GNU_SOURCE \
	-DHAVE_ISNAN=1 \
	-I/opt/local/include \
	-I/usr/include/gdal \
	-Iauxs

LIBS= \
	-L/opt/local/lib \
	-L/usr/lib

WARNS=-Wno-missing-braces \
      -Wno-unused-variable \

DEBUG=$(WARNS) -Wall -ggdb -O0
#DEBUG=$(WARNS) -O2

DEPEND_CPP=$(patsubst %.cpp,%.o,$(wildcard auxs/*.cpp)) $(patsubst %.cpp,%.o,$(wildcard *.cpp)) 
DEPEND_C_MAIN=$(patsubst %.c,%.o,$(wildcard auxs/*.c)) $(patsubst %.c,%.o,$(wildcard *.c))
DEPEND_C=$(filter-out main.o join.o, $(DEPEND_C_MAIN))

all: main join

main: $(DEPEND_CPP) $(DEPEND_C) main.o
	g++ -std=c++11 $(DEBUG) $(INCLUDES) $(DEPEND_C) $(DEPEND_CPP) main.o -o main $(LIBS) -lgdal -lgeos -lgeos_c

join: $(DEPEND_CPP) $(DEPEND_C) join.o
	g++ -std=c++11 $(DEBUG) $(INCLUDES) $(DEPEND_C) $(DEPEND_CPP) join.o -o join $(LIBS) -lgdal -lgeos -lgeos_c

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
