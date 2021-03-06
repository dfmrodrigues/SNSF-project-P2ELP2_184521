clc
close all
clear
profile on

k1 = 0.053;
k2 = 0.128;
k3 = 0.028;
cBin = 5.0;
cA0 = 0.72;
cB0 = 0.05;
cC0 = 0.08;
cD0 = 0.01;
V0 = 1.0;
cBmax = 0.025;
cDmax = 0.15;
tfmax = 250.0;
Fmin = 0.0;
Fmax = 2.0;
T = [cA0*V0,-1,0,0,0;cB0*V0,-1,-2,-1,cBin;cC0*V0,1,0,0,0;cD0*V0,0,1,0,0;V0,0,0,0,1];
nA0 = T(1,1);
nB0 = T(2,1);
nC0 = T(3,1);
nD0 = T(4,1);
V0 = T(5,1);
nAx = num2cell(T(1,2:end));
nBx = num2cell(T(2,2:end));
nCx = num2cell(T(3,2:end));
nDx = num2cell(T(4,2:end));
Vx = num2cell(T(5,2:end));
glhs = str2sym({'nA','nB','nC','nD','V'});
grhs = str2sym({'nA0+nAx{1}*x1+nAx{2}*x2+nAx{3}*x3+nAx{4}*x4',...
    'nB0+nBx{1}*x1+nBx{2}*x2+nBx{3}*x3+nBx{4}*x4',...
    'nC0+nCx{1}*x1+nCx{2}*x2+nCx{3}*x3+nCx{4}*x4',...
    'nD0+nDx{1}*x1+nDx{2}*x2+nDx{3}*x3+nDx{4}*x4',...
    'V0+Vx{1}*x1+Vx{2}*x2+Vx{3}*x3+Vx{4}*x4'});
phi_form = subs(str2sym({'-nC';'nB-cBmax*V';'nD-cDmax*V';'((t-tfmax)/50)^5+((t-tfmax)/50)^3+(t-tfmax)/50'}),glhs,grhs);
f_form = subs(str2sym({'k1*nA*nB/V';'k2*nB^2/V';'k3*nB';'u1/1e3';'u2'}),glhs,grhs);
x0 = str2sym({'0';'0';'0';'0'});
lb = str2sym({'Fmin'});
ub = str2sym({'Fmax'});
gh_form = str2sym({});
h_form = str2sym({});
clhs = str2sym({'u1';'u2'});
crhs = {subs(str2sym({'x5';'1e3*(k1^2*nA*nB^2*(-2*nA*nB+nB^2+6*cBin*nA*V-2*cBin*nB*V)+2*cBin*V*(8*k2^2*nB^4+4*cBin^2*(x5/1000)^2*V^2-cBin*nB*(x5/1000)*V*(5*(x5/1000)+7*k3*V)+2*k2*nB^2*(4*nB*(x5/1000)+5*k3*nB*V-6*cBin*(x5/1000)*V)+nB^2*((x5/1000)^2+5*k3*(x5/1000)*V+3*k3^2*V^2))-2*k1*nA*nB*(2*k2*nB^2*(nB-5*cBin*V)+V*(k3*nB*(nB-6*cBin*V)+cBin*(x5/1000)*(-5*nB+7*cBin*V))))/(2*cBin*nB*V^2*(-nB+cBin*V))'}),glhs,grhs)};
nci = [];
np = 1;
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

np = 0;
ni = 1;
idx_x0 = length(x0)+(1:ni);
av = [-2,1,-1];
ti = zeros(1,0);
tsc = [5.57,230.01];
tfc = 250;
x0c = [x0;1.299];
p0c = -1.1279e-3;
sc = [10,10,10,1,1e-3];
hessian_check(nci,1:nphi,av,idx_x0,ti,tsc,tfc,x0c,p0c,sc,deriv);

tf = 250;
av = [-2,1,-1];
nt = 2;%3;
scv = [10,10,10,1,1e-3];
deriv = 0;
n = 6;
m = 1000;
cheby = false;
suppl = 0;
xunits = 'min';
yunits = 'mL min$^{-1}$';
[tau,p_sol,d_sol,deg,tel,fvalv,uv,N,p_mon,p_s,d,dus,flag,cphi,As,ys,y,Av,yv,u0,du,sc] = hessian_solve(tf,x0,nci,nphi,av,np,ni,nt,ti,scv,deriv,lb,ub,n,m,cheby,suppl,xunits,yunits);

profile viewer