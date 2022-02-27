
INCLUDES= \
	-D_GNU_SOURCE \
	-DHAVE_ISNAN=1 \
	-I/opt/local/include \
	-I/usr/include/gdal \
	-Iinclude \
	-Irtree

LIBS= \
	-L/opt/local/lib \
	-L/usr/lib

WARNS=-Wno-missing-braces \
      -Wno-unused-variable \

DEBUG=$(WARNS) -Wall -ggdb -O0
#DEBUG=$(WARNS) -O2

SRC_C=$(patsubst %.c,%.o,$(wildcard rtree/*.c))
SRC_CPP=$(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
EXPERIMENTS=$(patsubst %.cpp,%.exp,$(wildcard experiments/*.cpp))

all: $(EXPERIMENTS)

%.exp: $(SRC_C) $(SRC_CPP) %.o
	clang++ -std=c++17 $(DEBUG) $(INCLUDES) $^ -o $@ $(LIBS) -lgdal -lgeos -lgeos_c
	ln -sf $@ $*

%.o: %.c
	clang -std=c11 $(DEBUG) $(INCLUDES) -c $< -o $@

%.o: %.cpp
	clang++ -std=c++17 $(DEBUG) $(INCLUDES) -c $< -o $@

clean:
	rm -f src/*.o
	rm -f rtree/*.o
	rm -f experiments/*.o
	rm -f experiments/*.exp

.PRECIOUS: %.cpp %.o

