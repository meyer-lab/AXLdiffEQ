clc; clear;

load 3wider1MH58R

names ={'xRev1','xFwd3','xRev3','AXLint2',...
    'kDeg','AXL2','Gas1'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

cutoff = min(fitt)+3;

params(fitt > cutoff,:) = [];
fitt(fitt > cutoff) = [];

[fitt, IDX] = sort(fitt);
params = params(IDX,:);


%%
clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

figure(1);

for ii = 1:7
    subplot(3,3,ii);
    
    semilogy(params(:,ii), fitt,'.','MarkerSize',10);
    title(names(ii));
    axis([minn(ii) maxx(ii) min(fitt) max(fitt)]);
end

