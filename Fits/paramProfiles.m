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

params(fitt > min(fitt)+1,:) = [];
fitt(fitt > min(fitt)+1) = [];

clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

params = [ones(size(params,1),1)*0.06, 10.^params];
params(:,end) = params(:,end) > 0.05;

%%


eOpt = @(x) ensembleOpt(x,params);

opts = psoptimset('Display','iter','TimeLimit',60*5,'CompletePoll','on');

vecc = [0 0 3 0 3 0 0 0 0 0 0 0 0 0];


pp = patternsearch(eOpt,zeros(14,1),[],[],[],[],-vecc,vecc,[], opts);
pp = pp';

%%

for ii = 1:size(params,1)
    pY = cLib_profile (linspace(0, 10, 100), params(ii,:).*(10.^pp), 10, 0);

    if ~isempty(pY)
        pY = pY / pY(1);
        plot(linspace(0, 10, 100), pY);
        axis([0 10 0 5]);
        hold on;
    end
end

