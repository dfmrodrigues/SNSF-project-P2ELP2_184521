function [Q0,z,tau,status,si] = sdp_mosek(s0_A,s0_barA,s0_b,n_p,q,w,u0k_s,uck_s,Q0_s,Qc_s,p,n_c)
[subk0,subl0,subj0] = ind2sub([u0k_s,u0k_s,n_p*p],1:(Q0_s*n_p));
[subkc,sublc,subjc] = ind2sub([uck_s,uck_s,n_p*n_c],1:(Qc_s*n_p));
flagsub0 = (subk0>=subl0);
flagsubc = (subkc>=sublc);
subk = [subk0,subkc];
subl = [subl0,sublc];
subj = [subj0,n_p*p+subjc];
flagsub = [flagsub0,flagsubc];
subk = subk(flagsub);
subl = subl(flagsub);
subj = subj(flagsub);
if(q==0)
    si.SDP.c = zeros(1,0);
    cval = s0_barA(1,flagsub)';
    si.SDP.a = sparse(zeros(size(s0_A,1)-1,0));
    aval = s0_barA(2:end,flagsub)';
    si.SDP.blc = full(s0_b(2:end)');
    si.SDP.buc = full(s0_b(2:end)');
else
    si.SDP.c = -w;
    cval = zeros(nnz(flagsub),1);
    si.SDP.a = sparse(s0_A);
    aval = s0_barA(:,flagsub)';
    si.SDP.blc = full(s0_b(:)');
    si.SDP.buc = full(s0_b(:)');
end
si.SDP.bardim = [u0k_s*ones(1,n_p*p),uck_s*ones(1,n_p*n_c)];
[indjc,~,valc] = find(cval);
si.SDP.barc.subj = subj(indjc);
si.SDP.barc.subk = subk(indjc);
si.SDP.barc.subl = subl(indjc);
si.SDP.barc.val = valc';
[indja,india,vala] = find(aval);
si.SDP.bara.subi = india';
si.SDP.bara.subj = subj(indja);
si.SDP.bara.subk = subk(indja);
si.SDP.bara.subl = subl(indja);
si.SDP.bara.val = vala';
si.param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-9;
si.param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-9;
si.param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-9;
[~,si.res] = mosekopt('minimize info',si.SDP,si.param);
if(q==0)
    z = w*[-1;si.res.sol.itr.y];
    tau = -cval'*si.res.sol.itr.barx+s0_b(1);
else
    z = si.res.sol.itr.y;
    tau = si.res.sol.itr.xx(1:q);
end
barQ0_s = u0k_s*(u0k_s+1)/2*p;
barQ0 = reshape(si.res.sol.itr.barx(1:barQ0_s*n_p),u0k_s*(u0k_s+1)/2,1,n_p*p);
ind0 = sub2ind([u0k_s,u0k_s,n_p*p],subk0(flagsub0),subl0(flagsub0),subj0(flagsub0));
Q0 = zeros(u0k_s,u0k_s,n_p*p);
Q0(ind0) = barQ0;
for m = 1:n_p
    for k = 1:p
        Q0(:,:,(m-1)*p+k) = Q0(:,:,(m-1)*p+k)+tril(Q0(:,:,(m-1)*p+k),-1)';
    end
end
status = si.res.sol.itr.solsta;
end