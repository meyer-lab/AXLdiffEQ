function fit_full_sensitivity_analysisKD(ii)

params = log10(getOptimParams(1));

names = {'B1','B2','U1','U2','xFwd1','xRev1','xFwd3','xRev3',...
    'AXLint1','AXLint2','scaleA','kRec','kDeg','fElse','fD2',...
    'Gas1','Gas2','AXL1','AXL2'};

minn = log10([1.2,1E-6,0.042,1E-5,1E-5,1E-5,1E-5,1E-5,1E-6,1E-3,1E-6,...
    1E-3,1E-4,1E-5,1E-1,1E-6,1E-6,1,1]);

maxx = log10([1.2,0.6,0.042,1E5,1E5,1E5,1E5,1E5,1,100,1,0.1,0.1,1,1,...
    100,100,1E5,1E5]);

N = 25;

A = zeros(1,19);

if ii == 1
    A(4) = 1; A(2) = -1;
elseif ii == 2
    A(6) = 1; A(5) = -1;
elseif ii == 3
    A(8) = 1; A(7) = -1;
end

curr = A*params';
vv = linspace(curr-2,curr+2,N);

Fopts = optimoptions('fmincon','Algorithm','interior-point','MaxFunEvals',1e6);
Dopts = psoptimset('Vectorized','on','TimeLimit',1600);

for jj = 1:N
    [~,fitIDXlocal(jj)] = fmincon(@cLib,params,[],[],A,vv(jj),minn,maxx,[],Fopts);
    [~,fitIDXglobal(jj)] = patternsearch(@cLib,params,[],[],A,vv(jj),minn,maxx,[],Dopts);
 
    save(['KD' mat2str(ii)]);
end