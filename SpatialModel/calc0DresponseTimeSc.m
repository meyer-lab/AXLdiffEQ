function ZZ = calc0DresponseTimeSc(params)

TPS = [0 logspace(-2,4,30) 1e6];

detail = 5;

Gass = logspace(-4,0.0969,detail);
AXXL = logspace(1,4,detail);
autoC = logspace(-9,-1,detail); % Up to -1

for ii = 1:length(AXXL)
    for jj = 1:length(Gass)    
        for xx = 1:length(autoC)
            if ii == 1 && jj == 1 && xx == 1
                YY(1,:) = cLib_profile (TPS, params, AXXL(ii), autoC(xx), Gass(jj));
                condition(1,:) = [AXXL(ii) autoC(xx) Gass(jj)];
            else
                YY(end+1,:) = cLib_profile (TPS, params, AXXL(ii), autoC(xx), Gass(jj));
                condition(end+1,:) = [AXXL(ii) autoC(xx) Gass(jj)];
            end
            
            YY(end,:) = abs(YY(end,:) - YY(end,end));
            
            YY(end,:) = YY(end,:) / YY(end,1);
        end
    end
    
end

ZZ = median(YY);