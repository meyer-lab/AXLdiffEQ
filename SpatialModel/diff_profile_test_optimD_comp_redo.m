% Diff profile test
% This script generated Figure 2D
clc; clear;

names = {'Gas','A','A1','A2','A12','D1','D2','Ai','A1i','A2i','A12i','D1i','D2i','Gasi'};
ttt = [1 0 0 0 0 0 0 0 0 0 0 0 0 1;...
       0 1 0 0 0 0 0 1 0 0 0 0 0 0;...
       0 0 1 0 0 0 0 0 1 0 0 0 0 0;...
       0 0 0 1 0 0 0 0 0 1 0 0 0 0;...
       0 0 0 0 1 0 0 0 0 0 1 0 0 0;...
       0 0 0 0 0 1 0 0 0 0 0 1 0 0;...
       0 0 0 0 0 0 1 0 0 0 0 0 1 0];

xx = linspace(0,1,90);
sshape = @(x) cos(xx*pi/3).^x(1);

Gass = @(x) 2*sshape(x)/mean(sshape(x).*xx)*(10^x(2));
meanify = @(x) 2*x'*xx'/length(xx);

D = ones(1,14)*10;
D([1 3:end]) = 0;

ffSpat = @(x,y) cLib_diff_profile (30, 359.46, 0.0581, Gass(x), D, 1, 1);

A = logspace(-2,3,20);

for ii = 1:length(A)
    temp = meanify(ffSpat([A(ii) -1.204],1));
    B(ii,:) = ttt*temp;
    disp(ii);
end

for ii = 1:size(B,2)
    C(:,ii) = B(:,ii) - B(1,ii);
end

semilogx(A,C);
legend(names);