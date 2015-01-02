function outter = cLibA549bothJustB (in)

in2 = 10.^in;

if in(13) > 0.5
    in2(13) = 1;
else
    in2(13) = 0;
end

in7 = [0.06 in2(1:13)];
in4 = in7;
in4(2) = in2(14);
in4(4) = in2(15);

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

outter = calllib('libOptimize','pyEntry',libpointer('doublePtr',in7));

outter = outter + calllib('libOptimize','pyEntryFull',libpointer('doublePtr',in4));