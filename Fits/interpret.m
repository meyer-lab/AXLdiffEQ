clc; clear;

load widertEXKsu
fitStruct2 = fitStruct;
load widerY4h4NgrH
fitStruct = [fitStruct, fitStruct2];

names ={'U2', 'xFwd1','xRev1','xFwd3','xRev3', 'AXLint1','AXLint2',...
    'kRec','kDeg','fElse','AXL2','Gas1','picker'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

min(fitt(params(:,end) > 0.05))

min(fitt(params(:,end) < 0.05))


cutoff = min(fitt)+2;

params(fitt > cutoff,:) = [];
fitt(fitt > cutoff) = [];

IDXrem = params(:,7) < params(:,6);

%fitt(IDXrem) = [];
%params(IDXrem,:) = [];

[fitt, IDX] = sort(fitt);
params = params(IDX,:);











%%
clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

%figure(1);

% subplot(2,3,1);
% scatter(params(:,2),params(:,3),20,fitt,'filled');
% xlabel('xFwd_1');
% ylabel('xRev_1');
% axis([minn(2) maxx(2) minn(3) maxx(3)]);
% colormap('parula')
% brighten(-0.6)
% 
% subplot(2,3,2);
% scatter(params(:,4),params(:,5),20,fitt,'filled');
% xlabel('xFwd_3');
% ylabel('xRev_3');
% axis([minn(4) maxx(4) minn(5) maxx(5)]);
% 
% subplot(2,3,3);
% scatter(params(:,6),params(:,7),20,fitt,'filled');
% xlabel('AXLint1');
% ylabel('AXLint2');
% axis([minn(6) maxx(6) minn(7) maxx(7)]);
% 
% subplot(2,3,4);
% scatter(params(:,8),params(:,9),20,fitt,'filled');
% xlabel('kRec');
% ylabel('kDeg');
% axis([minn(8) maxx(8) minn(9) maxx(9)]);
% 
% subplot(2,3,5);
% scatter(params(:,11),params(:,12),20,fitt,'filled');
% xlabel('AXL2');
% ylabel('Gas1');
% axis([minn(11) maxx(11) minn(12) maxx(12)]);
% 
% subplot(2,3,6);
% scatter(params(:,1),params(:,10),20,fitt,'filled');
% xlabel('U2');
% ylabel('fElse');
% axis([minn(1) maxx(1) minn(10) maxx(10)]);




%%

params = [ones(size(params,1),1)*0.06, 10.^params];
params(:,end) = params(:,end) > 0.05;
% 
% xx = linspace(0,1,80);
% GasConc = 100;
% 
% GassF = @(x, au) (xx < x)*GasConc + au/100; 
% 
% D = zeros(1,14);
% D(2) = 1000;
% D(7:end) = 0.1;
% 
% tps = linspace(0,10,100);
% 
% %figure(2);
% 
% for ii = 1:size(params,1)
%     B = cLib_diff_profile_pYavg (tps, params(ii,:), GassF(0.10, params(ii,end-1)), D, 1);
% 
%     try
%         outter(ii) = B(end) / B(1); %#ok<SAGROW>
%         
%         plot(tps, B / B(1));
%         hold on;
%         drawnow;
%     catch
%         outter(ii) = NaN;
%     end
%     
% end

%%

% figure(2);
% 
% Gass = logspace(log10(0.25),log10(64),20);
% 
% for ii = 1:size(params,1)
%     ii
%     
%     for jj = 1:length(Gass)
%         stimProfile = cLib_profile ([0 60 240], params(ii,:), Gass(jj), 0);
%         
%         stimProfileTot = cLib_profile ([0 60 240], params(ii,:), Gass(jj), 2);
%         
%         stimProfileSurf = cLib_profile ([0 60 240], params(ii,:), Gass(jj), 3);
% 
%         if numel(stimProfile) == 0
%             outHour(ii,jj) = nan; %#ok<*SAGROW>
%             outFour(ii,jj) = nan;
%             
%             outHourSurf(ii,jj) = nan;
%             outFourSurf(ii,jj) = nan;
%             
%             outHourTot(ii,jj) = nan;
%             outFourTot(ii,jj) = nan;
%             continue;
%         end
% 
%         outHour(ii,jj) = stimProfile(2) / stimProfile(1);
%         outFour(ii,jj) = stimProfile(3) / stimProfile(1);
%         
%         outHourTot(ii,jj) = stimProfileTot(2) / stimProfileTot(1);
%         outFourTot(ii,jj) = stimProfileTot(3) / stimProfileTot(1);
%         
%         outHourSurf(ii,jj) = stimProfileSurf(2) / stimProfileSurf(1);
%         outFourSurf(ii,jj) = stimProfileSurf(3) / stimProfileSurf(1);
%     end
% end
% 
% subplot(1,3,1);
% semilogx(Gass,outHour','b');
% hold on;
% semilogx(Gass,outFour','k');
% axis([min(Gass) max(Gass) 0 max(1,max(max([outHour outFour])))]);
% 
% subplot(1,3,2);
% semilogx(Gass,outHourTot','b');
% hold on;
% semilogx(Gass,outFourTot','k');
% axis([min(Gass) max(Gass) 0 max(1,max(max([outHourTot outFourTot])))]);
% 
% subplot(1,3,3);
% semilogx(Gass,outHourSurf','b');
% hold on;
% semilogx(Gass,outFourSurf','k');
% axis([min(Gass) max(Gass) 0 max(1,max(max([outHourSurf outFourSurf])))]);
% 


%%


%plot(log10(params(:,2) ./ params(:,1)),fitt,'.')

Kd = (params(1,2) ./ params(1,1));
Kd1 = 0.035;
x = logspace(2, 5,1000);

xData = [64 16 4 1 0.25];
yData = [3443.11 3143.41 3018.88 2608.88 2690.24];


plotyy(log10(x/1000), Kd ./ (Kd + x),log10(xData), yData)








