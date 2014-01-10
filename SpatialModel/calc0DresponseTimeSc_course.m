clc; clear;

params = getOptimParams(1);

TPS = linspace(0,10,1000);

detail = 20;

Gass = logspace(-2,0.0969,detail);
AXXL = logspace(-1,3,detail);
autoC = logspace(-4,0,detail);

for ii = 1:length(AXXL)
    for jj = 1:length(Gass)    
        for xx = 1:length(autoC)
            if ii == 1 && jj == 1 && xx == 1
                YY(1,:) = cLib_profile (TPS, params, AXXL(ii), autoC(xx), Gass(jj), 0);
                condition(1,:) = [AXXL(ii) autoC(xx) Gass(jj)];
            else
                YY(end+1,:) = cLib_profile (TPS, params, AXXL(ii), autoC(xx), Gass(jj), 0);
                condition(end+1,:) = [AXXL(ii) autoC(xx) Gass(jj)];
            end
            
            YY(end,:) = YY(end,:) / YY(end,1);
        end
    end
    
end

figure(1);
plot(TPS,YY')





%% Long time scale
clc; clear;
TPS = 240;

clear YY;

Gass = logspace(-2,1,200);
AXXL = logspace(-1,3,5);
autoC = 0.001; % Up to -1

color = 'kbrg';

for xx = 1
    
    params = getOptimParams(xx);
    
    for ii = 1:length(AXXL)
        for jj = 1:length(Gass)
                YY(jj) = cLib_profile (TPS, params, AXXL(ii), autoC, Gass(jj), 1);
        end
        
        YY = YY / YY(1);
        
        if ii > 1 || xx > 1
            hold on;
        end
        semilogx(Gass, YY',color(ii));
        drawnow;
        
    end
end

