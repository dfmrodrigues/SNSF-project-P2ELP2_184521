clc
close all
clear
profile on

tfmax = 12.0;
umin = 0.0;
umax = 1.0;
glhs = str2sym({'s1','s2','q'});
grhs = str2sym({['x',num2str(1)],['x',num2str(2)],['x',num2str(3)]});
phi_form = [subs(str2sym({'q'}),glhs,grhs);str2sym('((t-tfmax)/2)^5+((t-tfmax)/2)^3+((t-tfmax)/2)^1');str2sym('((tfmax-t)/2)^5+((tfmax-t)/2)^3+((tfmax-t)/2)^1')];
f_form = subs(str2sym({'s1-s1*s2-0.4*s1*u1';'-s2+s1*s2-0.2*s2*u1';'(s1-1)^2+(s2-1)^2+1e-4*u1^2'}),glhs,grhs);
x0 = str2sym({'0.5';'0.7';'0.0'});
lb = str2sym({'umin'});
ub = str2sym({'umax'});
gh_form = str2sym({});
h_form = str2sym({});
clhs = str2sym({'u1'});
crhs = cell(1,0);
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
nphi = size(phi_form,1);

idx_x0 = length(x0)+(1:ni);
av = [-1,-2,-1];
ti = zeros(1,0);
tsc = [3,6];
tfc = 12;
x0c = x0;
p0c = zeros(0,1);
sc = [0.3,1,1];
hessian_check(nci,1:nphi,av,idx_x0,ti,tsc,tfc,x0c,p0c,sc,deriv);

tf = 12;
av = [-1,-2,-1];
nt = 2;%3;
scv = [0.3,1,1];
deriv = 5e-1/3;
n = 10;
m = 2000;
cheby = false;
suppl = 0;
[tau,p_sol,d_sol,deg,tel,fvalv,uv,N,p_mon,p_s,d,dus,flag,cphi,As,ys,y,Av,yv,u0,du,sc] = hessian_solve(tf,x0,nci,nphi,av,np,ni,nt,ti,scv,deriv,lb,ub,n,m,cheby,suppl);

profile viewer