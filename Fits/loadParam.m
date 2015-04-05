function [fitt, params, names, params2, cutoff, minn, maxx] = loadParam ()

load WithpYnewBalance_FP7;
% fitStruct2 = fitStruct; %#ok<NODEF>
% load WithpYnewBalance_56v
% fitStruct = [fitStruct, fitStruct2];
% fitStruct2 = fitStruct;
% load WithpYnewBalance_B2a;
% fitStruct = [fitStruct, fitStruct2];
% fitStruct2 = fitStruct;
% load WithpYnewBalance_3YD;
% fitStruct = [fitStruct, fitStruct2];

names = {'U2','xFwd1','xRev4','kDeg','fPhase','AXL','Gas','pD'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

disp(min(fitt(params(:,end) > 0.5)));
disp(min(fitt(params(:,end) < 0.5)));

cutoff = min(fitt)+3;

params(fitt > cutoff,:) = [];
fitt(fitt > cutoff) = [];

[fitt, IDX] = sort(fitt);
params = params(IDX,:);
params = flipud(params);
fitt = fliplr(fitt);

params2 = 10.^params;
params2(:,end) = params2(:,end) > 0.05;