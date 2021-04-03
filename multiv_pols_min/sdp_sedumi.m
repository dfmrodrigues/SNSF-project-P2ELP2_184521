function [Q0,z,tau,status,si] = sdp_sedumi(s0_A,s0_barA,s0_b,n_p,q,w,u0k_s,uck_s,Q0_s,Qc_s,p,n_c)
if(q==0)
    si.SDP.A = s0_barA(2:end,:);
    si.SDP.b = s0_b(2:end);
    si.SDP.c = s0_barA(1,:)';
else
    si.SDP.A = [s0_A,s0_barA];
    si.SDP.b = s0_b(:);
    si.SDP.c = sparse([-w';zeros((Q0_s+Qc_s)*n_p,1)]);
end
si.SDP.K.f = q;
si.SDP.K.s = [u0k_s*ones(1,n_p*p),uck_s*ones(1,n_p*n_c)];
si.param.eps = 1e-9;
si.param.bigeps = eps^(1/4);
si.param.errors = 1;
[si.res.x,si.res.y,si.res.info] = sedumi(si.SDP.A,si.SDP.b,si.SDP.c,si.SDP.K,si.param);
if(q==0)
    z = w*[-1;si.res.y];
    tau = -si.SDP.c'*si.res.x+s0_b(1);
else
    z = si.res.y;
    tau = si.res.x(1:q);
end
Q0 = reshape(si.res.x((q+1):(q+Q0_s*n_p)),u0k_s,u0k_s,n_p*p);
if(si.res.info.numerr==2)
    status = "Failed";
elseif(si.res.info.err(end)>eps^(1/2))
    status = "Inaccurate/Solved";
else
    status = "Solved";
end
end