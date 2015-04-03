function [pY, tot, surf] = cLib_profile (tps, params, GasStim)

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

out1 = libpointer('doublePtr',1:length(tps));
out2 = libpointer('doublePtr',1:length(tps));
out3 = libpointer('doublePtr',1:length(tps));
tpsIn = libpointer('doublePtr',tps);
paramsIn = libpointer('doublePtr',params);

x = calllib('libOptimize','calcProfileMatlab',out1,out2,out3,paramsIn,tpsIn,uint32(length(tps)),GasStim);

if x == 0
    pY = out1.Value;
    tot = out2.Value;
    surf = out3.Value;
else
    pY = [];
    tot = [];
    surf = [];
end
