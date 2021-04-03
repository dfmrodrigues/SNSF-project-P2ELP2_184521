clc
close all
clear
profile on

m1 = 38*2.39;
m2 = 0.003*0.82;
Tinf = 298.15; % K
Tb = 308.15; % K
Rho = 2800.0; % kg m-3
Cp = 795.0; % J kg-1 K-1
d = 0.2e-3; % m
r = 1.5e-3; % m
eta = 0.4;
K = 1.43; % W m-2 K-1
T0 = 310.15; % K
Pmin = 1.0; % W
Pmax = 5.0; % W
Tmax = 316.15; % K
CEMsp = 1.5; % min
k1 = m2*eta/Rho/Cp/pi/r^2/d;
k2 = 2*pi*r*d*K*m1/Rho/Cp/pi/r^2/d*(1/2+(Tmax-Tb)/(Tb-Tinf))*log(1+(Tb-Tinf)/(Tmax-Tb));

glhs = str2sym({'CEM','T','P'});
grhs = str2sym({'x1','x2','u1'});
phi_form = subs(str2sym({'t/100';'T-T0';'CEMsp-CEM'}),glhs,grhs);
f_form = subs(str2sym({'0.5^(43.0+273.15-T)/60.0';'k1*P-k2*((T-Tinf)-(T-Tb))/log((T-Tinf)/(T-Tb))'}),glhs,grhs);
gh_form = str2sym({});
h_form = subs(str2sym({'T-Tmax'}),glhs,grhs);
x0 = str2sym({'0.0';'T0'});
lb = str2sym({'Pmin'});
ub = str2sym({'Pmax'});
clhs = str2sym({'u1'});
crhs = {subs(str2sym({'k2*((T-Tinf)-(T-Tb))/log((T-Tinf)/(T-Tb))/k1'}),glhs,grhs)};
nci = 1;
np = 0;
ni = 0;
nt = 2;
nti = 0;
phi_form = vpa(subs(phi_form));
f_form = vpa(subs(f_form));
gh_form = vpa(subs(gh_form));
h_form = vpa(subs(h_form));
x0 = double(subs(x0));
lb = double(subs(lb));
ub = double(subs(ub));
clhs = vpa(subs(clhs));
for i = 1:length(crhs)
    crhs{i} = vpa(subs(crhs{i}));
end
deriv = 2;

[nci,phi,dphidxt,d2phidxt2,f,dfdxp,d2Hdxp2,d2fdxp2,h,dhdx] = hessian_compile(nci,clhs,crhs,lb,ub,phi_form,f_form,gh_form,h_form,np,ni,nt,nti,deriv);
nphi = length(phi_form);

idx_x0 = length(x0)+(1:ni);
av = [-2,-1];
ti = zeros(1,0);
tsc = 92.8404967453808;
tfc = 109.4495929298785;
x0c = x0;
p0c = [];
sc = [100,100];
hessian_check(nci,1:nphi,av,idx_x0,ti,tsc,tfc,x0c,p0c,sc,deriv);

tf = 120;
av = [-2,-1];
nt = 2;
scv = [100,100];
deriv = 1e-2;
n = 6;
m = 1000;
cheby = false;
suppl = 0;
[tau,p_sol,d_sol,deg,tel,fvalv,uv,N,p_mon,p_s,d,dus,flag,cphi,As,ys,y,Av,yv,u0,du,sc] = hessian_solve(tf,x0,nci,nphi,av,np,ni,nt,ti,scv,deriv,lb,ub,n,m,cheby,suppl);

profile viewer