clc; clear;

load WithpY_Xa8;
fitStruct2 = fitStruct;
load WithpY_nRN;
fitStruct = [fitStruct, fitStruct2];

names = {'U2','xFwd1','xRev3','AXLint1','AXLint2','kRec','kDeg','fElse','AXL','Gas'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

cutoff = min(fitt)+6;

params(fitt > cutoff,:) = [];
fitt(fitt > cutoff) = [];

[fitt, IDX] = sort(fitt);
params = params(IDX,:);

params = 10.^params;
params(:,end) = params(:,end) > 0.05;

%%
clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

tps = [240];
Gass = logspace(-2, log10(64), 10);

for ii = 1:size(params,1)
    ppp = params(ii,:);
    
    for jj = 1:length(Gass)
        for xx = 1:12
             temp(xx) = cLib_profile (tps, ppp, Gass(jj), xx-1);
        end
        
        stimProfile(ii,jj,1:6) = temp(1:6) + temp(7:12);
    end
    
    ii/size(params,1)
    
    loglog(Gass,stimProfile(ii,:,1) / stimProfile(ii,1,1),'r');
    hold on;
    loglog(Gass,stimProfile(ii,:,2) / stimProfile(ii,1,2),'k');
    drawnow;
end


%axis([min(tps) max(tps) 0 100]);