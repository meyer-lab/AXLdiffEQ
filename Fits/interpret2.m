clc; clear;

load 2wider4P73jF

names ={'U2','xRev1','xFwd3','xRev3', 'AXLint1','AXLint2',...
    'kDeg','fElse','AXL2','Gas1','picker'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

cutoff = min(fitt)+3;

params(fitt > cutoff,:) = [];
fitt(fitt > cutoff) = [];

IDXrem = params(:,7) < params(:,6);

fitt(IDXrem) = [];
params(IDXrem,:) = [];

[fitt, IDX] = sort(fitt);
params = params(IDX,:);


%%
clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

figure(1);

for ii = 1:10
    subplot(2,5,ii);
    
    semilogy(params(:,ii), fitt,'.','MarkerSize',10);
    title(names(ii));
    axis([minn(ii) maxx(ii) min(fitt) max(fitt)]);
end


