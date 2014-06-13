function outter = cLibA549endoVary (in)

in = 10.^in;

in = [in zeros(size(in,1),1)];

%unloadlibrary('libOptimize');
if ~libisloaded('libOptimize')
    loadlibrary('libOptimize.dylib','BlasHeader.h')
end

in2 = reshape(in',1,numel(in));

out = libpointer('doublePtr',1:size(in,1));
in3 = libpointer('doublePtr',in2);


x = calllib('libOptimize','matlabEntryA549VaryEndo',out,in3,size(in,1));

if x == 0
    outter = out.Value;
else
    outter = 1e6;
end