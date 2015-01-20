clc; clear;

load noPYMq0;
% fitStruct2 = fitStruct;
% load widenoPYxLI;
% fitStruct = [fitStruct, fitStruct2];
% load widenoPYe4M;
% fitStruct = [fitStruct, fitStruct2];

names = {'U2','xFwd1','xRev3','AXLint1','AXLint2','kRec','kDeg','fElse','AXL','Gas'};

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

params2 = 10.^params;
params2(:,end) = params(:,end) > 0.05;

%%
clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

figure(1);

subplot(2,3,1);
scatter(params(:,2),params(:,3),25,fitt,'filled');
xlabel(names(2));
ylabel(names(3));
axis([minn(2) maxx(2) minn(3) maxx(3)]);
colormap('gray')
brighten(-0.2)

subplot(2,3,2);
scatter(params(:,4),params(:,5),25,fitt,'filled');
xlabel(names(4));
ylabel(names(5));
axis([minn(4) maxx(4) minn(5) maxx(5)]);

subplot(2,3,3);
scatter(params(:,6),params(:,7),25,fitt,'filled');
xlabel(names(6));
ylabel(names(7));
axis([minn(6) maxx(6) minn(7) maxx(7)]);

subplot(2,3,4);
scatter(params(:,9),params(:,10),20,fitt,'filled');
xlabel(names(9));
ylabel(names(10));
axis([minn(9) maxx(9) minn(10) maxx(10)]);

subplot(2,3,5);
scatter(params(:,1),params(:,8),20,fitt,'filled');
xlabel(names(1));
ylabel(names(8));
axis([minn(1) maxx(1) minn(8) maxx(8)]);

%%

params = 10.^params;
params(:,end) = params(:,end) > 0.05;

figure(2);

Gass = [logspace(log10(0.25),log10(64),10)];

for ii = 1:size(params,1)
    for jj = 1:length(Gass)
        stimProfile = cLib_profile ([0 60 240], params(ii,:), Gass(jj), 0);
        
        stimProfileTot = cLib_profile ([0 60 240], params(ii,:), Gass(jj), 2);
        
        stimProfileSurf = cLib_profile ([0 60 240], params(ii,:), Gass(jj), 3);

        if numel(stimProfile) == 0
            outHour(ii,jj) = nan; %#ok<*SAGROW>
            outFour(ii,jj) = nan;
            
            outHourSurf(ii,jj) = nan;
            outFourSurf(ii,jj) = nan;
            
            outHourTot(ii,jj) = nan;
            outFourTot(ii,jj) = nan;
            continue;
        end

        outHour(ii,jj) = stimProfile(2) / stimProfile(1);
        outFour(ii,jj) = stimProfile(3) / stimProfile(1);
        
        outHourTot(ii,jj) = stimProfileTot(2) / stimProfileTot(1);
        outFourTot(ii,jj) = stimProfileTot(3) / stimProfileTot(1);
        
        outHourSurf(ii,jj) = stimProfileSurf(2) / stimProfileSurf(1);
        outFourSurf(ii,jj) = stimProfileSurf(3) / stimProfileSurf(1);
    end
end

subplot(1,3,1);
semilogx(Gass,outHour','b');
hold on;
semilogx(Gass,outFour','k');
axis([min(Gass) max(Gass) 0 1.5]);

subplot(1,3,2);
semilogx(Gass,outHourTot','b');
hold on;
semilogx(Gass,outFourTot','k');
axis([min(Gass) max(Gass) 0 1.4]);

subplot(1,3,3);
semilogx(Gass,outHourSurf','b');
hold on;
semilogx(Gass,outFourSurf','k');
axis([min(Gass) max(Gass) 0 1.5]);



figure(3);
Gass = 1.25;
tps = linspace(0,10,100);

for ii = 1:size(params,1)
        stimProfile2{ii} = cLib_profile (tps, params(ii,:), Gass, 0);
        
        plot(tps,stimProfile2{ii} / stimProfile2{ii}(1));
        hold on;
end





%%


% plot(log10(params(:,2) ./ params(:,1)),fitt,'.')
% 
% Kd = (params(1,2) ./ params(1,1));
% Kd1 = 0.035;
% x = logspace(2, 5,1000);
% 
% xData = [64 16 4 1 0.25];
% yData = [3443.11 3143.41 3018.88 2608.88 2690.24];
% 
% 
% plotyy(log10(x/1000), Kd ./ (Kd + x),log10(xData), yData)








