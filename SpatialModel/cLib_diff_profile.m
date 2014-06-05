function outter = cLib_diff_profile (tps, params, AXLin, autocrine, GasIn, Din, endoImpair, degImpair)

Nspecies = 14;

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

dataPtr = libpointer('doublePtr',1:(length(tps)*length(GasIn)*Nspecies));
GasInPtr = libpointer('doublePtr',GasIn);
pIn = libpointer('doublePtr',params);
pTps = libpointer('doublePtr',tps);
DinP = libpointer('doublePtr',Din);


x = calllib('libOptimize','matlabDiffTPS',dataPtr, AXLin, GasInPtr, length(GasIn), autocrine, pIn, pTps, length(tps), DinP, endoImpair, degImpair);

if x == 0
    outter = dataPtr.Value;
    
    outter = reshape(outter,length(GasIn),[],length(tps));
    
else
    disp('Error!');
    outter = -1*ones([length(GasIn) Nspecies length(tps)]);
end
