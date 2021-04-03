clc
close all
clear
profile on

phi_form = str2sym({'x1';'x2';'-x2';'x3';'-x3';'t-2'});
f_form = str2sym({'1/2*x2^2';'x3+u1';'-u1'});
x0 = str2sym({'0';'1';'1'});
lb = str2sym({'-10'});
ub = str2sym({'10'});
gh_form = str2sym({});
h_form = str2sym({});
clhs = str2sym({'u1'});
crhs = {str2sym({'-x2-x3'})};
nci = [];
np = 0;
ni = 0;
nt = 3;
nti = 0;
phi_form = vpa(subs(phi_form));
f_form = vpa(subs(f_form));
x0 = double(subs(x0));
lb = double(subs(lb));
ub = double(subs(ub));
gh_form = vpa(subs(gh_form));
h_form = vpa(subs(h_form));
clhs = vpa(subs(clhs));
for i = 1:length(crhs)
    crhs{i} = vpa(subs(crhs{i}));
end
deriv = 2;

[nci,phi,dphidxt,d2phidxt2,f,dfdxp,d2Hdxp2,d2fdxp2,h,dhdx] = hessian_compile(nci,clhs,crhs,lb,ub,phi_form,f_form,gh_form,h_form,np,ni,nt,nti,deriv);
nphi = length(phi_form);

idx_x0 = length(x0)+(1:ni);
av = [-2,1,-1];
ti = zeros(1,0);
tsc = [0.298965977,1.927246292];
tfc = 2.0;
x0c = x0;
p0c = [];
sc = [1,1,1];
hessian_check(nci,1:nphi,av,idx_x0,ti,tsc,tfc,x0c,p0c,sc,deriv);

lam = zeros(nphi,nt+length(x0)+1);
pi0 = zeros(nt+length(x0),1);
[u,fval,~,~,nu] = hessian_fmincon(nci,nphi,av,np,idx_x0,ti,tsc,tfc,x0,[],sc,lam,pi0,lb,ub,0.001);
ts = u(1:2);
tf = u(3);
addpath('./hessian');
hessian_plot(nci,1:nphi,av,ti,0,ts,tf,x0,[],lam,pi0,nu.ineqnonlin)
rmpath('./hessian');

profile viewer