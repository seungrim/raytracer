CC = g++
DEBUG = 
CFLAGS = -c -O3 $(DEBUG)
LFLAGS = $(DEBUG)

all:raytracer 

raytracer: raytracer.o rtraytracer.o rtmultilayer.o mlcalc.o rtsurface.o rtshape.o fminbnd.o
	$(CC) $(LFLAGS) raytracer.o rtraytracer.o rtmultilayer.o mlcalc.o rtsurface.o rtshape.o fminbnd.o -o raytracer

raytracer.o: raytracer.cpp rtraytracer.h rtsurface.h rtshape.h rtmultilayer.h rtnklibrary.h
	$(CC) $(CFLAGS) raytracer.cpp

rtraytracer.o: rtraytracer.cpp rtraytracer.h rtsurface.h
	$(CC) $(CFLAGS) rtraytracer.cpp

rtsurface.o: rtsurface.cpp rtsurface.h rtshape.h fminbnd.h rtmultilayer.h mlcalc.h
	$(CC) $(CFLAGS) rtsurface.cpp

rtshape.o: rtshape.cpp rtshape.h fminbnd.h
	$(CC) $(CFLAGS) rtshape.cpp

rtmultilayer.o: rtmultilayer.cpp rtmultilayer.h mlcalc.h
	$(CC) $(CFLAGS) rtmultilayer.cpp

mlcalc.o: mlcalc.cpp mlcalc.h
	$(CC) $(CFLAGS) mlcalc.cpp

fminbnd.o: fminbnd.cpp fminbnd.h
	$(CC) $(CFLAGS) fminbnd.cpp

clean:
	rm -rf *o raytracer 