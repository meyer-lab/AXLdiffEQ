function interpret()

load WithpYnewBalance_Z48;
fitStruct2 = fitStruct; %#ok<NODEF>
load WithpYnewBalance_56v
fitStruct = [fitStruct, fitStruct2];
fitStruct2 = fitStruct;
load WithpYnewBalance_B2a;
fitStruct = [fitStruct, fitStruct2];

names = {'U2','xFwd1','xRev4','int1','int2','kRec','kDeg','fElse','AXL','Gas'};

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

params = 10.^params;
params(:,end) = params(:,end) > 0.05;

figure(2);

Gass = [0 logspace(log10(0.25),log10(64),100)];

for ii = 1:size(params,1)
    
    for jj = 1:length(Gass)
            [stimProfile, stimProfileTot, stimProfileSurf] = cLib_profile ([0 60 240], params(ii,:), Gass(jj));

            outHour(ii,jj) = stimProfile(2);
            outFour(ii,jj) = stimProfile(3);

            outHourTot(ii,jj) = stimProfileTot(2);
            outFourTot(ii,jj) = stimProfileTot(3);

            outHourSurf(ii,jj) = stimProfileSurf(2);
            outFourSurf(ii,jj) = stimProfileSurf(3);
    end
end

Gass(1) = 0.1;

for ii = 1:size(params,1)

    subplot(1,3,1);
    semilogx(Gass,outHour(ii,:)','Color',[1 1 1] + [1 1 0]*(fitt(ii) - max(fitt))/(cutoff - min(fitt)));
    hold on;
    semilogx(Gass,outFour(ii,:)','Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)));
    axis([min(Gass) max(Gass) 0 max(outFour(ii,:))*1.1]);

    subplot(1,3,2);
    semilogx(Gass,outHourTot(ii,:)','Color',[1 1 1] + [1 1 0]*(fitt(ii) - max(fitt))/(cutoff - min(fitt)));
    hold on;
    semilogx(Gass,outFourTot(ii,:)','Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)));
    axis([min(Gass) max(Gass) 0 max(outFourTot(ii,:))*1.1]);

    subplot(1,3,3);
    semilogx(Gass,outHourSurf(ii,:)','Color',[1 1 1] + [1 1 0]*(fitt(ii) - max(fitt))/(cutoff - min(fitt)));
    hold on;
    semilogx(Gass,outFourSurf(ii,:)','Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)));
    axis([min(Gass) max(Gass) 0 max(outFourSurf(ii,:))*1.1]);

end



figure(3);
Gass = 1.25;
tps = linspace(0,10,100);

for ii = 1:size(params,1)
        stimProfile2 = cLib_profile (tps, params(ii,:), Gass);
        
        plot(tps,stimProfile2,'Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)));
        axis([min(tps) max(tps) 0 14]);
        hold on;
end

end


function plotData (names, minn, maxx, params, fitt, idx1, idx2)

scatter(params(:,idx1),params(:,idx2),25,fitt,'filled');
xlabel(names(idx1));
ylabel(names(idx2));
axis([minn(idx1) maxx(idx1) minn(idx2) maxx(idx2)]);
colormap('gray')
brighten(-0.2)

end

