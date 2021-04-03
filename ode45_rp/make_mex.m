function make_mex
eval('mex -output rp45_mex rp45_mex.c rp45.c ntrp45.c');
eval('mex -output ntrp45_mex ntrp45_mex.c ntrp45.c');
eval('mex -output rp3h_mex rp3h_mex.c rp3h.c ntrp3h.c');
eval('mex -output ntrp3h_mex ntrp3h_mex.c ntrp3h.c');
eval('mex -output binterp_mex binterp_mex.c binterp.c');
end
