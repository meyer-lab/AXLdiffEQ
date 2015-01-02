
clc; clear;


params = [-1.292 -0.5406 -2.0000 -2.9996 -1.2757 -0.2799 4.1348 -1.3263 0];
params = 10.^params;

parF = [0.06 6];

outter = [];

for ii = 1:5000
    pp = [parF 10.^(-60 + 60*rand) 10.^(0.799) 10.^(-60 + 60*rand) params];
    
    
    B = cLib_profile ([0 10], pp, 20, 1);
    
    
    if isempty(B)
        continue;
    end
    
    if B(2) / B(1) > 10
        outter(end+1,:) = log10([pp([3 5]), B(2) / B(1)]);
    end
end

[~, idx] = sort(outter(:,end));


outter = outter(idx,:);

%%
subplot(1,2,1);
plot(outter(:,1),outter(:,2),'o')
subplot(1,2,2);
plot(outter(:,1)-outter(:,2),'o')
