INCLUDES= \
	-I/opt/local/include \
	-Iauxs

LIBS=-L/opt/local/lib

all: geosext.o
	gcc -ggdb -O0 -c $(INCLUDES) main.c histogram.c auxs/*.c 
	g++ -ggdb -O0 *.o -o main $(LIBS) -lgdal -lgeos -lgeos_c


geosext.o: 
	g++ -ggdb -O0 -std=c++11 -c $(INCLUDES) auxs/*.cpp

clean:
	rm *.o
