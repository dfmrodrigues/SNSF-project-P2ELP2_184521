% An implementation of direct single shooting
% Adapted from Joel Andersson, 2016

function [J_opt,T_opt,x_opt] = calc_single_shooting(tf,N,x0_opt,u_opt,p,phi_opt,psi_opt,F_opt)

import casadi.*

x_opt = zeros(length(x0_opt),N+1);
x_opt(:,1) = x0_opt;
J_opt = 0;
for k = 1:N
    fk = F_opt{k}('x0',x_opt(:,k),'p',p,'u',u_opt(:,k));
    x_opt(:,k+1) = full(fk.xf);
    J_opt = J_opt+full(fk.qf);
end
phif = phi_opt('xf',x_opt(:,end),'tf',tf);
J_opt = J_opt+full(phif.phif);
psif = psi_opt('xf',x_opt(:,end),'tf',tf);
T_opt = full(psif.psif);

end
