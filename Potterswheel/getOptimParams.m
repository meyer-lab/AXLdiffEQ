function params3 = getOptimParams(idxSel)

if isempty(idxSel)
    idxSel = 1;
end

load new_cochran2_noBT;

aRR = [3 4 5 6 7 8 9 10 17 13 14 11 12 16 19 15 18] - 2;

% Clear scales
nameID = strfind(names,'obs1');
for ii = length(nameID):-1:1
    if ~isempty(nameID{ii})
        params(:,ii) = [];
        names(ii) = [];
    end
end

[~, idx] = sort(chiSq);

%nnn = nnn(idx);
%chiSq = chiSq(idx);
params = params(idx,aRR);
%names = names(aRR);


params3 = params(idxSel,:);

params3 = [1.2, params3(1), 0.042, params3(2:end)];