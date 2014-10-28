function outter = cLib_profile (params)

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

paramsIn = libpointer('doublePtr',params);

outter = calllib('libOptimize','calcLemke',paramsIn);

