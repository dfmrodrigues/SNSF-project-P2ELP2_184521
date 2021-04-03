function make_mex
eval('mex -output uniquerows_mex uniquerows_mex.c heapify.c uniquesorted.c');
eval('mex -output ismemberrows_mex ismemberrows_mex.c');
end
