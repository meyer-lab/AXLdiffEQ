function [pY, tot, surf, convF, species] = cLib_profile (tps, params, GasStim)

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

out1 = libpointer('doublePtr',1:length(tps));
out2 = libpointer('doublePtr',1:length(tps));
out3 = libpointer('doublePtr',1:length(tps));
out4 = libpointer('doublePtr',1:(length(tps)*13));
convFac = libpointer('doublePtr',1:3);

x = calllib('libOptimize','calcProfileMatlab',out1,out2,out3,out4,...
    libpointer('doublePtr',params),libpointer('doublePtr',tps),uint32(length(tps)),GasStim,convFac);

if x == 0
    pY = out1.Value;
    tot = out2.Value;
    surf = out3.Value;
    convF = convFac.Value;
    species = out4.Value;
else
    pY = [];
    tot = [];
    surf = [];
    convF = [];
    species = [];
end
