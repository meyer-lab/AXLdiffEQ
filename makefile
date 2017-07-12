.PHONY: all clean

DEPS = Code/libOptimize/BlasHeader.h Code/libOptimize/ModelRunning.h Code/libOptimize/CVodeHelpers.h Code/libOptimize/cobyla.h

%.o: %.c $(DEPS)
	$(CXX) -c -I. -std=c++0x -Wall -fPIC -o $@ $<

%.dylib: %.so
	cp $< $@

Code/libOptimize/cobyla.o: Code/libOptimize/cobyla.c
	gcc -c -o $@ $<

Code/libOptimize.so: Code/libOptimize/BlasHeader.o Code/libOptimize/ModelRunning.o Code/libOptimize/CVodeHelpers.o Code/libOptimize/cobyla.o
	cd Code/libOptimize; g++ -shared -fPIC -o ../libOptimize.so *.o -I. -lsundials_nvecserial -lsundials_cvode

Code/libOptimize.h: Code/libOptimize/BlasHeader.h
	grep "(double" $< > $@

all: Code/libOptimize.so Code/libOptimize.h Code/libOptimize.dylib

clean:
	rm -f Code/libOptimize/*.o Code/libOptimize.so Code/libOptimize.h Code/libOptimize.dylib