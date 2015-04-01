function interpret()

load WithpYnewBalance_Z48;

names = {'U2','xFwd1','xRev4','int1','int2','kRec','kDeg','fElse','AXL','Gas'};

fitt = [];
params = [];

for ii = 1:length(fitStruct)
    fitt = [fitt, fitStruct{ii}.fitIDXglobal];
    params = [params; fitStruct{ii}.paramOpt];
end

disp(min(fitt(params(:,end) > 0.5)));
disp(min(fitt(params(:,end) < 0.5)));

cutoff = min(fitt)+6;

params(fitt > cutoff,:) = [];
fitt(fitt > cutoff) = [];

[fitt, IDX] = sort(fitt);
params = params(IDX,:);
params = flipud(params);
fitt = fliplr(fitt);

params2 = 10.^params;
params2(:,end) = params(:,end) > 0.05;


clear fitIDXglobal fitStruct slices symbols A b Dopts ii xxxx fname paramOpt

figure(1);

plotD = @(idx1, idx2) plotData(names, minn, maxx, params, fitt, idx1, idx2);

subplot(2,3,1);
plotD(2, 3); % 

subplot(2,3,2);
plotD(1, 5);

subplot(2,3,3);
plotD(4, 6); %

subplot(2,3,4);
plotD(9, 10); %

subplot(2,3,5);
plotD(7, 8); %

%%

% params = 10.^params;
% params(:,end) = params(:,end) > 0.05;
% 
% figure(2);
% 
% Gass = logspace(log10(0.25),log10(64),100);
% 
% for ii = 1:size(params,1)
%     
%     surfFrac(ii) = cLib_profile (0, params(ii,:), 0, 4) / cLib_profile (0, params(ii,:), 0, 2) / 135.2;
%     
%     for jj = 1:length(Gass)
%             stimProfile = cLib_profile ([0 60 240], params(ii,:), Gass(jj), 1);
% 
%             stimProfileTot = cLib_profile ([0 60 240], params(ii,:), Gass(jj), 2);
% 
%             stimProfileSurf = cLib_profile ([0 60 240], params(ii,:), Gass(jj), 4);
% 
%             outHour(ii,jj) = stimProfile(2) / stimProfile(1);
%             outFour(ii,jj) = stimProfile(3) / stimProfile(1);
% 
%             outHourTot(ii,jj) = stimProfileTot(2) / stimProfileTot(1);
%             outFourTot(ii,jj) = stimProfileTot(3) / stimProfileTot(1);
% 
%             outHourSurf(ii,jj) = stimProfileSurf(2) / stimProfileSurf(1);
%             outFourSurf(ii,jj) = stimProfileSurf(3) / stimProfileSurf(1);
%     end
% end
% 
% 
% for ii = 1:size(params,1)
% 
%     subplot(1,3,1);
%     semilogx(Gass,outHour(ii,:)' / mean(outHour(ii,:)),'Color',[1 1 1] + [1 1 0]*(fitt(ii) - max(fitt))/(cutoff - min(fitt)));
%     hold on;
%     semilogx(Gass,outFour(ii,:)' / mean(outFour(ii,:)),'Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)));
%     axis([min(Gass) max(Gass) 0 3]);
% 
%     subplot(1,3,2);
%     semilogx(Gass,outHourTot(ii,:)','Color',[1 1 1] + [1 1 0]*(fitt(ii) - max(fitt))/(cutoff - min(fitt)));
%     hold on;
%     semilogx(Gass,outFourTot(ii,:)','Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)));
%     axis([min(Gass) max(Gass) 0 1.4]);
% 
%     subplot(1,3,3);
%     semilogx(Gass,outHourSurf(ii,:)','Color',[1 1 1] + [1 1 0]*(fitt(ii) - max(fitt))/(cutoff - min(fitt)));
%     hold on;
%     semilogx(Gass,outFourSurf(ii,:)','Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)));
%     axis([min(Gass) max(Gass) 0 1.5]);
% 
% end
% 
% 
% 
% figure(3);
% Gass = 1.25;
% tps = linspace(0,10,100);
% 
% for ii = 1:size(params,1)
%         stimProfile2{ii} = cLib_profile (tps, params(ii,:), Gass, 0);
%         
%         plot(tps,stimProfile2{ii} / mean(stimProfile2{ii}),'Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)));
%         hold on;
% end
% 
% figure(4);
% bar(surfFrac);

end


function plotData (names, minn, maxx, params, fitt, idx1, idx2)

scatter(params(:,idx1),params(:,idx2),25,fitt,'filled');
xlabel(names(idx1));
ylabel(names(idx2));
axis([minn(idx1) maxx(idx1) minn(idx2) maxx(idx2)]);
colormap('gray')
brighten(-0.2)

end

