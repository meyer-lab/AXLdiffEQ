function outter = cLib_diff_profile_pY (tps, AXLin, autocrine, GasIn, Din, endoImpair, degImpair, frac)

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

dataPtr = libpointer('doublePtr',1:(length(tps)*length(GasIn)));
GasInPtr = libpointer('doublePtr',GasIn);
pIn = libpointer('doublePtr',getOptimParams(1));
pTps = libpointer('doublePtr',tps);
DinP = libpointer('doublePtr',Din);


x = calllib('libOptimize','matlabDiffTPS_pY',dataPtr, AXLin, GasInPtr, length(GasIn), autocrine, pIn, pTps, length(tps), DinP, endoImpair, degImpair, frac);

if x == 0
    outter = dataPtr.Value;
else
    disp('Error');
    outter = -1;
end

% extern "C" int matlabDiffTPS(double *dataPtr, double AXLin, double *GasIn, int gridIn, double autocrine, double *params, double *tps, int nTps)