INCLUDES= \
	-D_GNU_SOURCE \
	-I/opt/local/include \
	-I/usr/include/gdal \
	-Iauxs

LIBS= \
	-L/opt/local/lib \
	-L/usr/lib

all: geosext.o
	clang -std=c11 -ggdb -O3 -c $(INCLUDES) *.c auxs/*.c 
	clang++ -std=c++11 -ggdb -O3 *.o -o main $(LIBS) -lgdal -lgeos -lgeos_c


geosext.o: 
	clang++ -ggdb -O3 -std=c++11 -c $(INCLUDES) auxs/*.cpp

clean:
	rm *.o
