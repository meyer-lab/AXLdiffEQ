function outter = cLib_profile (tps, params, AXL, autocrine, GasStim, frac)

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

out = libpointer('doublePtr',1:length(tps));
tpsIn = libpointer('doublePtr',tps);
paramsIn = libpointer('doublePtr',params);

x = calllib('libOptimize','calcProfileMatlab',out,paramsIn,tpsIn,length(tps),autocrine,AXL,GasStim, frac);

if x == 0
    outter = out.Value;
else
    outter = [];
    disp('Failed eval');
end
