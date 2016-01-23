INCLUDES= \
	-I/opt/local/include \
	-Iauxs

LIBS=-L/opt/local/lib

all: geosext.o
	gcc -c $(INCLUDES) main.c histogram.c auxs/*.c 
	g++ *.o -o main $(LIBS) -lgdal -lgeos -lgeos_c


geosext.o: 
	g++ -std=c++11 -c $(INCLUDES) auxs/*.cpp

