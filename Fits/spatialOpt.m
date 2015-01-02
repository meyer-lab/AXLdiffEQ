function outter = spatialOpt(D)


params = [0.06, 6, 0.017323, 12.722, 0.96976, 1.9591, 0.25047, 0.020395,...
    0.031742, 0.081154, 0.38841, 13954, 0.0099596, 1];

GassF = @(au) (linspace(0,1,50) < D(1))*D(2) + au; 

D(1:2) = [];

B = cLib_diff_profile_pYavg ([0 10], params, GassF(params(end-1)), D, 1);


try
    outter = B(end) / B(1);
catch
    outter = 0;
end