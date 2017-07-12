% Performs the fitting procedure.
function A549SsensFullWider ()
    slices = 50;
    symbols = ['a':'z' 'A':'Z' '0':'9'];
    rng('shuffle');
    fname = symbols(randi(numel(symbols),[1 3]));

    minn = log10([0.6,1E-20, 1E-5,0.003,0.03,1E-3,1E-2, 1E-2,100,1E-3,   1]);
    maxx = log10([6E3,    1,  1E5,  0.3,   1,   1,   1,    1,1E5,   1,  10]);
                 %U2 , xFwd,xRev4, int1,int2,kRec,kDeg,fElse,AXL, Gas, pD1

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

        fitStruct{xxxx}.paramOpt = paramOpt; %#ok<AGROW>
        fitStruct{xxxx}.fitIDXglobal = fitIDXglobal; %#ok<AGROW>

        save(['fit_' fname]);
    end
end

% Calls the C routine to calculate error from parameter sets.
function outter = cLibLoc (in, N)
    in2 = 10.^in;
    in2(:,end) = in(:,end) > 0.5;

    if ~libisloaded('libOptimize')
        loadlibrary('libOptimize.dylib');
    end

    outP = libpointer('doublePtr',zeros(1,numel(in)/N));

    calllib('libOptimize','pyEntryVec',libpointer('doublePtr',in2'),outP,numel(in)/N);

    outter = outP.Value';
end