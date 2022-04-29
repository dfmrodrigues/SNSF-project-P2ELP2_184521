clc
close all
clear
profile on

phi_form = str2sym({'-x1';'0.6-x3';'t-0.25'});
f_form = str2sym({'x2';'(u1-310*(x2*(1+tanh(500*(x2+0.03)))/2)^2*exp(500*(1-x1)))/x3-1/x1^2';'-u1/0.5*(1+tanh(50*(x3-0.3)))/2'});
x0 = str2sym({'1';'0';'1'});
lb = str2sym({'0.0'});
ub = str2sym({'3.5'});
gh_form = str2sym({});
h_form = str2sym({});
clhs = str2sym({'u1'});
crhs = cell(1,0);
nci = [];
np = 0;
ni = 1;%2;
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
tsc = [0.02,0.07];
tfc = 0.2;
x0c = [x0;2];
p0c = 20;
sc = [0.1,0.1,0.1,1,10];
hessian_check(nci,1:nphi,av,idx_x0,ti,tsc,tfc,x0c,p0c,sc,deriv);

tf = 0.25;
av = 3;%[-2,1,-1];
nt = 3;
scv = [0.1,0.1,0.1,1,10];
deriv = 0;
n = 6;
m = 1000;
cheby = false;
suppl = 0;
[tau,p_sol,d_sol,deg,tel,fvalv,uv,N,p_mon,p_s,d,dus,flag,cphi,As,ys,y,Av,yv,u0,du,sc] = hessian_solve(tf,x0,nci,nphi,av,np,ni,nt,ti,scv,deriv,lb,ub,n,m,cheby,suppl);

profile viewer