clc; clear;

load WithpYnewBalance_3YD;

names = {'U2','xFwd1','xRev3','AXLint1','AXLint2','kRec','kDeg','fElse','AXL','Gas'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal]; %#ok<AGROW>
    params = [params; fitStruct{ii}.paramOpt]; %#ok<AGROW>
end

disp(min(fitt(params(:,end) > 0.5)));
disp(min(fitt(params(:,end) < 0.5)));

cutoff = min(fitt)+6;

params(fitt > cutoff,:) = [];
fitt(fitt > cutoff) = [];

[fitt, IDX] = sort(fitt);
params = params(IDX,:);

params = 10.^params;
params(:,end) = params(:,end) > 0.05;

clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

for ii = 1:size(params,1)

    Kd = (params(ii,1) ./ 0.06)/1000;
    Kd1 = 0.035/1000;
    x = logspace(-4, 0,1000);


    semilogx(x, Kd ./ (Kd + x),'Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)))
    hold on;
    semilogx(x, Kd1 ./ (Kd1 + x),'r');
end

axis([min(x) max(x) 0 1]);
xlabel('Gas6');
ylabel('Fraction free');
