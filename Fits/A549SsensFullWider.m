function A549SsensFullWider ()

slices = 50;
symbols = ['a':'z' 'A':'Z' '0':'9'];
rng('shuffle');
fname = symbols(randi(numel(symbols),[1 3]));

minn = log10([0.6 ,1E-20, 1E-5,1E-3, 1E-10,100,1E-3, 1]);
maxx = log10([6E3 ,    1,  1E5,   1,  1E10,1E5,   1,  10]);
             %U2  ,xFwd ,xRev4,kDeg,fPhase,AXL,Gas
             
Dopts = psoptimset('TimeLimit',45*60,'Display','off','CompletePoll','on',...
    'CompleteSearch','on','Vectorized','on');

parpool(11);

for xxxx = 1:100
    parfor ii = 0:(slices*length(maxx) - 1)
        IDX = mod(ii,length(maxx))+1;
        vv = linspace(minn(IDX),maxx(IDX),slices+1); %#ok<PFBNS>

        vIDX = floor(ii/length(maxx)) + 1;

        minn2 = minn; minn2(IDX) = vv(vIDX);
        maxx2 = maxx; maxx2(IDX) = vv(vIDX+1); 

        params = minn2 + (rand(size(minn2)) .* (maxx2 - minn2));
        
        try
            [paramOpt(ii+1,:),fitIDXglobal(ii+1)] = ...
                patternsearch(@(x) cLibLoc(x, length(minn)),params,[],[],[],[],...
                minn2,maxx2,[],Dopts);
        catch
            paramOpt(ii+1,:) = params;
            fitIDXglobal(ii+1) = 1E6;
        end
    end
    
    disp(xxxx);
    
    fitStruct{xxxx}.paramOpt = paramOpt; %#ok<AGROW>
    fitStruct{xxxx}.fitIDXglobal = fitIDXglobal; %#ok<AGROW>

    save(['WithpYnewBalance_' fname]);
end

end

function outter = cLibLoc (in, N)

in2 = 10.^in;
in2(:,end) = in(:,end) > 0.5;

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

outP = libpointer('doublePtr',zeros(1,numel(in)/N));

calllib('libOptimize','pyEntryVec',libpointer('doublePtr',in2'),outP,numel(in)/N);

outter = outP.Value';

end