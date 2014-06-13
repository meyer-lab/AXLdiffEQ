function A549stochFit
clc;
rng('shuffle');

symbols = ['a':'z' 'A':'Z' '0':'9'];
randName = symbols (randi(numel(symbols),[1 8]));

params = 10.^[-3.204 0.77799 -2.3861 3.5071 -4.9999 -0.9382 0.86126	0.92224	-2.0199	-2.9999	-1.0002	-2.7892	-1.0254	-1.6667	-0.00012853	-1.4651	2.7388];

minn = log10([0.0006,0.0006,1E-5,1E-5,... % 'B1','B2','U1','U2'
    1E-5,1E-5,1E-5,1E-5,... % 'xFwd1','xRev1','xFwd3','xRev3'
    1E-6,1E-3,1E-6,... % 'AXLint1','AXLint2','scaleA'
    1E-3,1E-4,1E-5,1E-1,1E-6,1]); % 'kRec','kDeg','fElse','fD2','Gas1','AXL2'

maxx = log10([6,6,1E5,1E5,... % 'B1','B2','U1','U2'
    1E5,1E5,1E5,1E5,... % 'xFwd1','xRev1','xFwd3','xRev3'
    1,10,1,... % 'AXLint1','AXLint2','scaleA'
    0.1,0.1,1,1,0.1,1E5]); % 'kRec','kDeg','fElse','fD2','Gas1','AXL2'
% ScaleA upped from 0.1

A = zeros(2,17);
A(1,3) = 1;
A(1,1) = -1;
A(2,4) = -1;
A(2,2) = 1;
b = [1; 0];


opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter');
problem = createOptimProblem('fmincon','objective',...
 @cLibA549,'x0',params,'Aineq',A,'bineq',b,'lb',minn,'ub',maxx,'options',opts);
gs = GlobalSearch('NumStageOnePoints',2000,'NumTrialPoints',4000);
[x,f] = run(gs,problem)