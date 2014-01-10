clc; clear;

reoorder = [1:8 11 12 15 16 17 18 19];
load('/Users/aaron/Desktop/fitFull/noXLINK_allSig.mat');
[a, idxx] = sort(chiSq);
chiSq = chiSq(idxx);
params = params(idxx,reoorder);
names = names(reoorder);


chiSq = chiSq(1000:-1:1);
params = params(1000:-1:1,:);

TPS = logspace(-2,4,32);

for ii = 1:length(chiSq)
    
    params2 = params(ii,:);
    
    calc0DresponseTimeSc(params2)
    
    semilogx(TPS,calc0DresponseTimeSc(params2));
    hold on;
    disp(['ChiSq: ' mat2str(chiSq(ii)) ' on ' mat2str(ii)]);
    drawnow;
end