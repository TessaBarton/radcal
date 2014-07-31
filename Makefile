CC=g++
LD=g++

OBJS=GdalFileIO.o radcal.o
LIBS=-lgdal -lm
INCLUDEDIRS= -I/usr/local/Cellar/eigen/3.2.1/include/eigen3

CCFLAGS= -O3 -g -p -c -Wall -std=c++11
LDFLAGS=-g -p

all: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) $(LIBS) -o radcal

radcal.o: radcal.cpp
	$(CC) $(CCFLAGS) $(INCLUDEDIRS) radcal.cpp -o radcal.o

GdalFileIO.o: GdalFileIO.cpp
	$(CC) $(CCFLAGS) $(INCLUDEDIRS) GdalFileIO.cpp -o GdalFileIO.o

debug: radcal
	gdb radcal

test: all
	./radcal


clean:
	rm -f *.o radcal
