clc; clear;

data = load('out.csv');
data = data(1:10:end,:);


tps = [0 10];
GasStim = 64;


parfor ii = 1:size(data,1)
    
        params = 10.^data(ii,:);
        outter = cLib_profile (tps, params, GasStim, 1);
        
        %params(end) = 0;
        
        %outter2 = cLib_profile (tps, params, GasStim, 0);
        try
            resultt(ii) = outter(2) / outter(1) ;
        catch
            resultt(ii) = -1;
        end
    
end

semilogy(resultt);
legend('show');