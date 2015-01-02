clc; clear;

load 77j87EXi

names ={'B2', 'U2', 'xFwd1','xRev1','xFwd3','xRev3', 'AXLint1','AXLint2',...
    'kRec','kDeg','fElse','AXL2','Gas1','picker'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

params(fitt > min(fitt)+3,:) = [];
fitt(fitt > min(fitt)+3) = [];

params = [ones(size(params,1),1)*0.06, 10.^params];
params(:,end) = params(:,end) > 0.05;

clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

options = saoptimset('Display','iter');

for ii = 1:size(params,1)
    optimF = @(x) optimFun(params(ii,:),x);
    
    xx(1:2,ii) = simulannealbnd(optimF,[4 1],[0 0],[4 4],options);
    
    
    plot(xx);
    colorbar;
    drawnow;
end





