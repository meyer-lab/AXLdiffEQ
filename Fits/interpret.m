function interpret()

[fitt, params, names, params2, cutoff, minn, maxx] = loadParam;

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

params = params2;

Gass = [0 logspace(log10(0.25),log10(64),10)];

for ii = 1:size(params,1)
    
    for jj = 1:length(Gass)
            [stimProfile, stimProfileTot, stimProfileSurf, convF, species] = ...
                cLib_profile ([0 60 240], params(ii,:), Gass(jj));

            outHour(ii,jj) = stimProfile(2)*convF(1);
            outFour(ii,jj) = stimProfile(3)*convF(1);

            outHourTot(ii,jj) = stimProfileTot(2);
            outFourTot(ii,jj) = stimProfileTot(3);

            outHourSurf(ii,jj) = stimProfileSurf(2)*convF(2);
            outFourSurf(ii,jj) = stimProfileSurf(3)*convF(2);
            
            surfFrac(ii) = stimProfileSurf(1) / stimProfileTot(1);
            
            speciesList(ii,:) = species(1:13);
    end
end

figure(2);
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
        [stimProfile2, ~, ~, convF] = cLib_profile (tps, params(ii,:), Gass);
        
        plot(tps,stimProfile2*convF(3),'Color',ones(1,3)*(fitt(ii) - min(fitt))/(cutoff - min(fitt)));
        axis([min(tps) max(tps) 0 14]);
        hold on;
end

figure(4)
bar(surfFrac);

figure(5)
imagesc(speciesList)

end


function plotData (names, minn, maxx, params, fitt, idx1, idx2)

scatter(params(:,idx1),params(:,idx2),25,fitt,'filled');
xlabel(names(idx1));
ylabel(names(idx2));
axis([minn(idx1) maxx(idx1) minn(idx2) maxx(idx2)]);
colormap('gray')
brighten(-0.2)

end

