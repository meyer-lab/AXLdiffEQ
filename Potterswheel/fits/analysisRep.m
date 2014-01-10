
clc; clear;

files = dir('*.mat');
pwClear;
xx = 1;


pwLoadRepository(files(1).name);

for ii = 2:length(files)
    h = waitbar(ii/length(files));
    
    ii
    
    pwAppendFitHistoryOfRepository(files(ii).name);
    
    if (mod(ii,10) == 0)
        [~, fitGroupsIDs, fitGroupsNFits, fitGroupsChisq] = pwFitHistoryGetFitGroups;
        
        for jj = 1:length(fitGroupsIDs) 
            pwFitHistoryDeleteFits(mean(fitGroupsChisq{jj} > 150)*100, jj);
        end
    end
end

close(h);

