function fit_full_sensitivity_analysisKD(jj)

params = log10(getOptimParams(1))';

names = {'B1','B2','U1','U2','xFwd1','xRev1','xFwd3','xRev3',...
    'AXLint1','AXLint2','scaleA','kRec','kDeg','fElse','fD2',...
    'Gas1','Gas2','AXL1','AXL2'};

minn = log10([0.0006,0.0006,1E-5,1E-5,1E-5,1E-5,1E-5,1E-5,1E-6,1E-3,1E-6,1E-3,1E-4,1E-5,1E-1,1E-6,1E-6,1,1]);
maxx = log10([6,6,1E5,1E5,1E5,1E5,1E5,1E5,1,100,1,0.1,0.1,1,1,100,100,1E5,1E5]);

N = 100;

A = zeros(1,19);

A(3) = 1; A(1) = -1;
A(4) = -1; A(2) = 1;

curr = A*params';
vv = linspace(-4,1,N);

Fopts = optimoptions('fmincon','Algorithm','interior-point','MaxFunEvals',1e6);
Dopts = psoptimset('Vectorized','on','TimeLimit',12*60*60);


[~,fitIDXlocal(jj)] = fmincon(@cLib,params,[],[],A,vv(jj),minn,maxx,[],Fopts);
[~,fitIDXglobal(jj)] = patternsearch(@cLib,params,[],[],A,vv(jj),minn,maxx,[],Dopts);

save(['KD_diff' mat2str(jj)]);