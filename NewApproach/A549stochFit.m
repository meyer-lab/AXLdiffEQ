%function A549stochFit
clc;
rng('shuffle');

minn = log10([0.0006,1E-5,1E-5,1E-5,1E-5,1E-5,... % ,'B2','U2','xFwd1','xRev1','xFwd3','xRev3'
    1E-6,1E-3,... % 'AXLint1','AXLint2',
    1E-3,1E-4,1E-5,1E-6,1]); % 'kRec','kDeg','fElse',,'Gas1','AXL2'

maxx = log10([6,1E5,1E5,1E5,1E5,1E5,... % ,'B2','U2','xFwd1','xRev1','xFwd3','xRev3'
    1,10,... % 'AXLint1','AXLint2',
    0.1,0.1,1,0.1,1E5]); % 'kRec','kDeg','fElse',,'Gas1','AXL2'

options = gaoptimset('Display','iter','Vectorized','on','PopulationSize',400);

ga(@cLib,length(maxx),[],[],[],[],minn,maxx,[],options)