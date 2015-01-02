function outter = cLib (in)

in = 10.^in;

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end


in3 = libpointer('doublePtr',in);


outter = calllib('libOptimize','pyEntry',in3);
