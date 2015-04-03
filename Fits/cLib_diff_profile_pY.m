function outter = cLib_diff_profile_pY (tps, params, GasIn, Din, frac)

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

dataPtr = libpointer('doublePtr',1:(length(tps)*length(GasIn)));
GasInPtr = libpointer('doublePtr',GasIn);
pIn = libpointer('doublePtr',params);
pTps = libpointer('doublePtr',tps);
DinP = libpointer('doublePtr',Din);

x = calllib('libOptimize','matlabDiffTPS_pY',dataPtr, GasInPtr, uint32(length(GasIn)), pIn, pTps, uint32(length(tps)), DinP, frac);

%matlabDiffTPS_pY(doublePtr, doublePtr, uint32, doublePtr, doublePtr, uint32, doublePtr, int32)

if x == 0
    outter = dataPtr.Value;
else
    disp('Error');
    outter = -1;
end

% extern "C" int matlabDiffTPS(double *dataPtr, double AXLin, double *GasIn, int gridIn, double autocrine, double *params, double *tps, int nTps)