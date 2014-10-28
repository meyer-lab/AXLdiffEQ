function outter = foldActivate (params, autocrine, AXL)

% This is a new approach

tps = 10;

GasStim = [0 25 50 100 250 750];

for ii = 1:length(GasStim)

    outter(ii) = cLib_profile (tps, params, AXL, autocrine, GasStim(ii), 0);
end