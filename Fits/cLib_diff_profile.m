function outter = cLib_diff_profile (tps, params, GasIn, Din)

Nspecies = 13;

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

dataPtr = libpointer('doublePtr',1:(length(tps)*length(GasIn)*Nspecies));
GasInPtr = libpointer('doublePtr',GasIn);
pIn = libpointer('doublePtr',params);
pTps = libpointer('doublePtr',tps);
DinP = libpointer('doublePtr',Din);


x = calllib('libOptimize','matlabDiffTPS',dataPtr, GasInPtr, length(GasIn), pIn, pTps, length(tps), DinP);

if x == 0
    outter = dataPtr.Value;
    
    outter = reshape(outter,length(GasIn),[],length(tps));
    
else
    disp('Error!');
    outter = -1*ones([length(GasIn) Nspecies length(tps)]);
end
