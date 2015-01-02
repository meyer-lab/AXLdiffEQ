function outter = cLibA549 (in)

in2 = 10.^in;

if in(end) > 0.5
    in2(end) = 1;
else
    in2(end) = 0;
end

in2 = [0.06 in2];

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

in3 = libpointer('doublePtr',in2);

outter = calllib('libOptimize','pyEntry',in3);