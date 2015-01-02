clc; clear;

load 77j87EXi

[~, IDX] = sort(params(:,15));
goodness = params(IDX,15);
params = params(IDX,1:14);

params(goodness > (min(goodness) + 1),:) = [];
goodness(goodness > (min(goodness) + 1),:) = [];




%%
sshape = @(x,xx) cos(xx*pi/3).^x(1);
xx = linspace(0,1,50);
GasConc = 10;

GassF = @(x, au) (xx < x)*GasConc + au; 

D = zeros(1,14);
D(2) = 1;
D(8:end) = 0.01;

%parfor_progress(size(params,1));

parfor ii = 1:size(params,1)
    
    B = cLib_diff_profile_pYavg ([0 10], 10.^params(ii,:), GassF(0.04, 10.^params(ii,end)), D, 1, 1, 0);
    %parfor_progress;
    try
        outter(ii) = B(2) / B(1); %#ok<SAGROW>
        
        B(2) / B(1)
    catch
        outter(ii) = NaN;
    end
    
    
end

%parfor_progress(0);