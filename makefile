.PHONY: all clean

CXX=g++ -fPIC
CXXFLAGS = -D_GLIBCXX_USE_CXX11_ABI=0 -I. -std=c++0x -fPIC
DEPS = Code/libOptimize/BlasHeader.h Code/libOptimize/ModelRunning.h Code/libOptimize/CVodeHelpers.h Code/libOptimize/cobyla.h

%.o: %.c $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Code/libOptimize/cobyla.o: Code/libOptimize/cobyla.c
	gcc -fPIC -c -o $@ $<

Code/libOptimize.so: Code/libOptimize/BlasHeader.o Code/libOptimize/ModelRunning.o Code/libOptimize/CVodeHelpers.o Code/libOptimize/cobyla.o
	cd Code/libOptimize; $(CXX) -shared -o ../libOptimize.so *.o -I. -lsundials_nvecserial -lsundials_cvode

Code/libOptimize.h: Code/libOptimize/BlasHeader.h
	grep "(double" $< > $@

all: Code/libOptimize.so Code/libOptimize.h

clean:
	rm -f Code/libOptimize/*.o Code/libOptimize.so Code/libOptimize.h
