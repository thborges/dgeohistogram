INCLUDES= \
	-I/usr/include \
	-I/usr/include/gdal \
	-Iauxs

LIBS=-L/usr/lib

all: geosext.o
	gcc -ggdb -std=c11  -O0 -c $(INCLUDES) main.c histogram.c auxs/*.c 
	g++ -ggdb -O0 -std=c++11 *.o -o main $(LIBS) -lgdal -lgeos -lgeos_c


geosext.o: 
	g++ -ggdb -O0 -std=c++11 -c $(INCLUDES) auxs/*.cpp

clean:
	rm *.o
