clc; clear;

load bothY4h4NgrH


names ={'U2', 'xFwd1','xRev1','xFwd3','xRev3', 'AXLint1','AXLint2',...
    'kRec','kDeg','fElse','AXL2','Gas1','picker','U2b','xFwd1b','xRev1b',...
    'xFwd3b','xRev3b'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

params(fitt > min(fitt)+4,:) = [];
fitt(fitt > min(fitt)+4) = [];

%params = [ones(size(params,1),1)*0.06, 10.^params];
%params(:,end) = params(:,end) > 0.05;

%%
clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

for ii = 1:18
    subplot(3,6,ii);
    
    plot(params(:,ii), fitt,'o');
    title(names(ii));
    axis([minn(ii) maxx(ii) min(fitt) max(fitt)]);
end


ppSet = [params(:,1:5), zeros(size(params,1),1); params(:,14:end), ones(size(params,1),1)];


ppSet = [ppSet(:,1:5), ppSet(:,3)-ppSet(:,2), ppSet(:,5)-ppSet(:,4), ppSet(:,6)];

%%


%SVMStruct = svmtrain(ppSet(:,[2 7]),logical(ppSet(:,8)), 'showplot',true);