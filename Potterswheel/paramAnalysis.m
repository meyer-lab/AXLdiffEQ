clc; clear;

load('out');


% Clear scales
nameID = strfind(names,'obs1');
for ii = length(nameID):-1:1
    if ~isempty(nameID{ii})
        params(:,ii) = [];
        names(ii) = [];
    end
end

[~, idx] = sort(chiSq);

nnn = nnn(idx);
chiSq = chiSq(idx);
params = params(idx,:);

idx = chiSq < 10;
nnn(idx) = [];
chiSq(idx) = [];
params(idx,:) = [];

min(chiSq)

clear idx

IDDD = 1:100;

Kd2 = params(:,2)./params(:,1);


subplot(2,2,1);
plot(log10(params(IDDD,1:6)));
legend(names(1:6));

subplot(2,2,2);

plot(log10(params(IDDD,7:10)));
legend(names(7:10));


subplot(2,2,3);
plot(log10(params(IDDD,11:end)));
legend(names(11:end));


subplot(2,2,4);
plot(chiSq(IDDD));

%plot(chiSq(1:300),log10(params(1:300,16)./params(1:300,15)),'o')


% % idxx = max(Kd1,Kd2) < 10;
% % goodID = find(idxx == 0,1,'first')
% 
% for ii = find(Kd1 > Kd2)
%     temp = Kd1(ii);
%     Kd1(ii) = Kd2(ii);
%     Kd2(ii) = temp;
%     
%     params(ii,14:20) = params(ii,[14 16 15 17 18 20 19]);
% end
% 
% 
% 
% 
% scatter(log10(Kd1(1:300)),log10(Kd2(1:300)),5,chiSq(1:300));

% params(goodID,:)'