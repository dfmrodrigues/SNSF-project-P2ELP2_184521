clc
close all
clear
profile on

k3 = 0.475783*exp(5.34-3760/323)*0.5/0.2;
k4 = 0.1;
Mw = diag([71,88.12,122.52,36.45,156.97,46])/1000;
Rho = diag([1093,859.17,1085.53,1486.26,1070,790])/1000;
T = [0,-1,-2,1000/71;13000,-1,-1,0;0,1,0,0;0,1,2,0;0,0,1,0;100000,0,0,0];
S = size(T,1);
Vmax = 10000;
umin = 0;
umax = 5;
nCl20 = T(1,1);
nBA0 = T(2,1);
nMBA0 = T(3,1);
V0 = ones(1,S)/Rho*Mw*T(:,1);
conv = 0.9*nBA0;
nCl2x = num2cell(T(1,2:end));
nBAx = num2cell(T(2,2:end));
nMBAx = num2cell(T(3,2:end));
Vx = num2cell(ones(1,S)/Rho*Mw*T(:,2:end));
glhs = str2sym({'nCl2','nBA','nMBA','V'});
grhs = str2sym({'nCl20+nCl2x{1}*x1+nCl2x{2}*x2+nCl2x{3}*x3',...
    'nBA0+nBAx{1}*x1+nBAx{2}*x2+nBAx{3}*x3',...
    'nMBA0+nMBAx{1}*x1+nMBAx{2}*x2+nMBAx{3}*x3',...
    'V0+Vx{1}*x1+Vx{2}*x2+Vx{3}*x3'});
phi_form = subs(str2sym({'t/1000';'(conv-nMBA)/1000'}),glhs,grhs);
f_form = subs(str2sym({'k3*nBA*nCl2/V';'k4*nCl2/V*k3*nBA*nCl2/V';'u1'}),glhs,grhs);
x0 = str2sym({'0';'0';'0'});
lb = str2sym({'umin'});
ub = str2sym({'umax'});
gh_form = str2sym({});
h_form = subs(str2sym({'V-Vmax'}),glhs,grhs);
clhs = str2sym({'u1'});
crhs = {subs(str2sym({'-(Vx{1}*(k3*nBA*nCl2/V)+Vx{2}*(k4*nCl2/V*k3*nBA*nCl2/V))/Vx{3}'}),glhs,grhs)};
nci = 1;
np = 0;
ni = 1;
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
av = [-2,1,-2];
ti = zeros(1,0);
tsc = [70,870];
tfc = 1307;
x0c = [x0;1.3];
p0c = 5e-04;
sc = [10,10,10,1,1e-3];
hessian_check(nci,1:nphi,av,idx_x0,ti,tsc,tfc,x0c,p0c,sc,deriv);

tf = 1500;
av = 3;%[-2,1,-2];
nt = 3;
scv = [10,10,10,1,1e-3];
deriv = 0;%1e-2;
n = 6;
m = 3000;
cheby = false;
suppl = 0;
[tau,p_sol,d_sol,deg,tel,fvalv,uv,N,p_mon,p_s,d,dus,flag,cphi,As,ys,y,Av,yv,u0,du,sc] = hessian_solve(tf,x0,nci,nphi,av,np,ni,nt,ti,scv,deriv,lb,ub,n,m,cheby,suppl);

profile viewer