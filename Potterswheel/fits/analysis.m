
clc; clear;

files = dir('*.mat');
xx = 1;

for ii = 1:length(files)
    load(files(ii).name);
    
    h = waitbar(ii/length(files));
    
    for jj = 1:length(pwGlobals2.fitHistory)
        params(xx,:) = pwGlobals2.fitHistory(jj).parsForFitEnd;
        chiSq(xx) = pwGlobals2.fitHistory(jj).optimization.chisq(end);
        nnn{xx} = files(ii).name;
        xx = xx + 1;
    end
    
    min(chiSq)
end

[a idx] = min(chiSq);
a

nnn{idx}

names = pwGlobals2.fitHistory(1).parsForFitIDsAndNames;

close(h);

save('out.mat','params','chiSq','names','nnn');