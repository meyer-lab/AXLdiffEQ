function outter = cLib_flux (species)

params = getOptimParams();



%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

out = libpointer('doublePtr',1:10);
specIn = libpointer('doublePtr',species);
paramsIn = libpointer('doublePtr',params);

x = calllib('libOptimize','matlabFluxify',out,specIn,paramsIn);

if x == 0
    outter = out.Value;
else
    outter = [];
    disp('Failed eval');
end
