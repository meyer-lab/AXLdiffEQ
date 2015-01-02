function outter = cLib_diff_profile_pYavg (tps, params, GasIn, Din, frac)


%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

dataPtr = libpointer('doublePtr',1:length(tps));
GasInPtr = libpointer('doublePtr',GasIn);

pIn = libpointer('doublePtr',params);
pTps = libpointer('doublePtr',tps);
DinP = libpointer('doublePtr',Din);


x = calllib('libOptimize','matlabDiffTPS_pYavg',dataPtr, GasInPtr, length(GasIn), pIn, pTps, length(tps), DinP, 1, 1, frac);

if x == 0
    outter = dataPtr.Value;
else
    outter = -1;
    disp('Error');
end



end