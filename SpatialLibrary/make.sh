#! /bin/bash


cd ./CVode/

for fname in *.c; do
	g++ -c -flto -fPIC -lm -march=core2 -O3 -o ${fname%%.*}.o $fname
done

cd ../

for fname in *.cpp; do
	g++ -c -flto -fPIC -lm -std=c++11 -O3 -I./CVode/ -o ${fname%%.*}.o -march=core2 $fname
done

for fname in *.c; do
        g++ -c -flto -fPIC -lm -std=c++11 -O3 -o ${fname%%.*}.o -march=core2 $fname
done

g++ -fPIC -Wl,--gc-sections -flto *.o ./CVode/*.o -shared -lstdc++ -o BlasHeader.so 

