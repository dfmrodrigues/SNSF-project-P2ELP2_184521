clc
close all
clear
profile on

tfmax = 12.0;
umin = 0.0;
umax = 1.0;
k = 0.2;
K = num2cell(k^2*eye(2));
betad = true;
fast = true;
chance = true;
if(~betad)
if(~fast)
w = num2cell([
    0.0201051368113866*ones(30,1)
    0.0119139564145941*ones(30,1)
    0.0012962636171673*ones(30,1)
    0.0000179764901854*ones(30,1)]);% Normally distributed prior, 120 points
theta = num2cell([
    0.160635659707143*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]';...
    0.373712306584429*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]';...
    0.602436406397929*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]';...
    0.866951943872305*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]']);% Normally distributed prior, 120 points
else
w = {
   0.210491910111102 0.210491910111102 0.210491910111102 0.210491910111102 0.027777777777778 0.027777777777778...
   0.027777777777778 0.027777777777778 0.011730312111120 0.011730312111120 0.011730312111120 0.011730312111120}';% Normally distributed prior, 12 points
theta = {
  -0.151386791613424  -0.151386791613424
  -0.151386791613424   0.151386791613424
   0.151386791613424  -0.151386791613424
   0.151386791613424   0.151386791613424
  -0.489897948556635   0.000000000000000
   0.489897948556635   0.000000000000000
   0.000000000000000  -0.489897948556635
   0.000000000000000   0.489897948556635
  -0.396335765891742  -0.396335765891742
  -0.396335765891742   0.396335765891742
   0.396335765891742  -0.396335765891742
   0.396335765891742   0.396335765891742};% Normally distributed prior, 12 points
end
else
if(~fast)
max_order = 15;
a = [2;2];
b = [2;2];
c = [1/2;1/2];
orth_pols = @(k,n)sqrt((2*n+a(k)+b(k)-1)*gamma(n+a(k)+b(k)-1)*factorial(n)/gamma(n+a(k))/gamma(n+b(k))*beta(a(k),b(k)))*jacobiP(n,a(k)-1,b(k)-1,sym('x')/c(k));
np = (max_order+1)/2;
thf = @(k)roots(double(coeffs(orth_pols(k,np),sym('x'),'All')));
th1gv = thf(1);
th2gv = thf(2);
[th1,th2] = meshgrid(th1gv,th2gv);
theta = num2cell([th1(:),th2(:)]);
wf = @(k,thetak)(2*np+a(k)+b(k)-1)./((c(k)^2-thetak.^2).*double(subs(diff(orth_pols(k,np),sym('x')),sym('x'),thetak)).^2);
w1gv = wf(1,th1gv);
w2gv = wf(2,th2gv);
[w1,w2] = meshgrid(w1gv,w2gv);
w = num2cell(w1(:).*w2(:));
else
w = {
   0.191044776119403
   0.106177150595849
   0.106177150595849
   0.106177150595849
   0.106177150595849
   0.047980869656593
   0.047980869656593
   0.047980869656593
   0.047980869656593
   0.030697498463597
   0.030697498463597
   0.016543986635707
   0.016543986635707
   0.016543986635707
   0.016543986635707
   0.016188049850202
   0.016188049850202
   0.016188049850202
   0.016188049850202};% Beta distributed prior, 19 points
theta = {
   0.000000000000000   0.000000000000000
  -0.164852899906490  -0.241242938928766
  -0.164852899906490   0.241242938928766
   0.164852899906490  -0.241242938928766
   0.164852899906490   0.241242938928766
  -0.310520825354915  -0.047549861429124
  -0.310520825354915   0.047549861429124
   0.310520825354915  -0.047549861429124
   0.310520825354915   0.047549861429124
   0.000000000000000   0.423168047213609
   0.000000000000000  -0.423168047213609
  -0.334619897649469  -0.407296116730353
  -0.334619897649469   0.407296116730353
   0.334619897649469  -0.407296116730353
   0.334619897649469   0.407296116730353
  -0.440109933805545  -0.210778895104562
  -0.440109933805545   0.210778895104562
   0.440109933805545  -0.210778895104562
   0.440109933805545   0.210778895104562};% Beta distributed prior, 19 points
end
end
sigma_z = num2cell(0.1);
m_z = 10000;
Z = num2cell(norminv((1/2/m_z):(1/m_z):(1-1/2/m_z),0,sigma_z{1})');
n_theta = size(theta,2);
m_theta = size(theta,1);
% glhs = str2sym({'s1','s2','s11','s21','s12','s22','th1','th2','s1l1','s2l1','s11l1','s21l1','s12l1','s22l1'});
% grhs = str2sym({'xl1_1','xl1_2','xl1_3','xl1_4','xl1_5','xl1_6','theta1l1','theta2l1','xl1_1','xl1_2','xl1_3','xl1_4','xl1_5','xl1_6'});
% phi_form = subs([str2sym({'-log(det(inv([K{1,1},K{1,2};K{2,1},K{2,2}])))-2*ph0_1'}),...
%     str2sym({'wl1/m_z*ph0_2'}),...
%     str2sym({['-(Zl2).''*inv(sigma_z{1}^2)*(Zl2)/2-log((2*pi)^(1/2)*det(sigma_z{1}^2)^(1/2))',...
%             '+(Zl2+[s21l1,s22l1]*[th1;th2]).''*inv(sigma_z{1}^2+[s21l1,s22l1]*[K{1,1},K{1,2};K{2,1},K{2,2}]*[s21l1,s22l1].'')*(Zl2+[s21l1,s22l1]*[th1;th2])/2',...
%             '+log((2*pi)^(1/2)*det(sigma_z{1}^2+[s21l1,s22l1]*[K{1,1},K{1,2};K{2,1},K{2,2}]*[s21l1,s22l1].'')^(1/2))']});...
%     str2sym('((t-tfmax)/2)^5+((t-tfmax)/2)^3+((t-tfmax)/2)^1'),str2sym('0'),str2sym('0');...
%     str2sym('((tfmax-t)/2)^5+((tfmax-t)/2)^3+((tfmax-t)/2)^1'),str2sym('0'),str2sym('0')],glhs,grhs);
% f_form = subs(str2sym({'s1-(1+th1)*s1*s2-0.4*s1*u1';'-s2+(1+th2)*s1*s2-0.2*s2*u1';...
%             '-s1*s2+(1-(1+th1)*s2-0.4*u1)*s11+(-(1+th1)*s1)*s21';'((1+th2)*s2)*s11+(-1+(1+th2)*s1-0.2*u1)*s21';...
%             '(1-(1+th1)*s2-0.4*u1)*s12+(-(1+th1)*s1)*s22';'s1*s2+((1+th2)*s2)*s12+(-1+(1+th2)*s1-0.2*u1)*s22'}),glhs,grhs);
% x0 = str2sym(repmat({'0.5';'0.7';'0';'0';'0';'0'},m_theta,1));
% lb = str2sym({'umin'});
% ub = str2sym({'umax'});
% gh_form = str2sym({});
% glhs = str2sym({'s1','s2','th1','th2','s1l1','s2l1','s1l2','s2l2','s1l3','s2l3'});
% grhs = str2sym({'xl1_1','xl1_2','theta1l1','theta2l1','xl1_1','xl1_2','xl2_1','xl2_2','xl3_1','xl3_2'});
% if(~betad)
%     phi0_form = str2sym({'-log(det(inv([K{1,1},K{1,2};K{2,1},K{2,2}])))-2*ph0_1'});
% else
%     phi0_form = str2sym({'-ph0_1'});
% end
% phi_form = subs([phi0_form,...
%     str2sym({'wl1/m_z*ph0_2'}),...
%     str2sym({'-(Zl2).''*inv(sigma_z{1}^2)*(Zl2)/2-log((2*pi)^(1/2)*det(sigma_z{1}^2)^(1/2))-log(ph0_3)'}),...
%     str2sym({'wl3*exp(-(Zl2+s2l1-s2l3).''*inv(sigma_z{1}^2)*(Zl2+s2l1-s2l3)/2)/((2*pi)^(1/2)*det(sigma_z{1}^2)^(1/2))'});...
%     str2sym('((t-tfmax)/2)^5+((t-tfmax)/2)^3+((t-tfmax)/2)^1'),str2sym('0'),str2sym('0'),str2sym('0');...
%     str2sym('((tfmax-t)/2)^5+((tfmax-t)/2)^3+((tfmax-t)/2)^1'),str2sym('0'),str2sym('0'),str2sym('0')],glhs,grhs);
% f_form = subs(str2sym({'s1-(1+th1)*s1*s2-0.4*s1*u1';'-s2+(1+th2)*s1*s2-0.2*s2*u1'}),glhs,grhs);
% x0 = str2sym(repmat({'0.5';'0.7'},m_theta,1));
% lb = str2sym({'umin'});
% ub = str2sym({'umax'});
% gh_form = str2sym({});
glhs = str2sym({'s1','s2','th1','th2','s1l1','s2l1','s1l2','s2l2','s1l3','s2l3'});
grhs = str2sym({'xl1_1','xl1_2','theta1l1','theta2l1','xl1_1','xl1_2','xl2_1','xl2_2','xl3_1','xl3_2'});
if(~betad)
    phi0_form = str2sym({'-log(det(inv([K{1,1},K{1,2};K{2,1},K{2,2}])))-2*ph0_1'});
else
    phi0_form = str2sym({'-ph0_1'});
end
phi_form = subs([phi0_form,...
    str2sym({'wl1/m_z*ph0_2'}),...
    str2sym({'-(Zl2).''*inv(sigma_z{1}^2)*(Zl2)/2-log((2*pi)^(1/2)*det(sigma_z{1}^2)^(1/2))-log(ph0_3)'}),...
    str2sym({['wl3*exp(-(((1+s2l1/10)*Zl2+s2l1-s2l3)/(1+s2l3/10)).''*inv(sigma_z{1}^2)*(((1+s2l1/10)*Zl2+s2l1-s2l3)/(1+s2l3/10))/2)',...
    '/(2*pi)^(1/2)/det(sigma_z{1}^2)^(1/2)/(1+s2l3/10)*(1+s2l1/10)']});...
    str2sym('((t-tfmax)/2)^5+((t-tfmax)/2)^3+((t-tfmax)/2)^1'),str2sym('0'),str2sym('0'),str2sym('0');...
    str2sym('((tfmax-t)/2)^5+((tfmax-t)/2)^3+((tfmax-t)/2)^1'),str2sym('0'),str2sym('0'),str2sym('0')],glhs,grhs);
f_form = subs(str2sym({'s1-(1+th1)*s1*s2-0.4*s1*u1';'-s2+(1+th2)*s1*s2-0.2*s2*u1'}),glhs,grhs);
x0 = str2sym(repmat({'0.5';'0.7'},m_theta,1));
lb = str2sym({'umin'});
ub = str2sym({'umax'});
gh_form = str2sym({});
if(~chance)
h_form = str2sym({});
clhs = str2sym({'u1'});
crhs = cell(1,0);
nci = [];
else
r = norminv(0.95);
xub = 6.6;
nf = length(f_form);
hl_form = '(s1-sub)';
dhldt0_form = '(s1-(1+th1)*s1*s2)';
dhldt1_form = '(0.4*s1)';
h_form{1} = '-(';
for l = 1:length(w)
    h_form{1} = [h_form{1},'+',num2str(w{l},15)];
    h_form{1} = [h_form{1},'*',replace(hl_form,{'s1','sub'},{['x',num2str(1+(l-1)*nf)],num2str(xub,15)})];
end
h_form{1} = [h_form{1},')^2+',num2str(r^2,15),'/',num2str(1+r^2,15),'*('];
for l = 1:length(w)
    h_form{1} = [h_form{1},'+',num2str(w{l},15)];
    h_form{1} = [h_form{1},'*',replace(hl_form,{'s1','sub'},{['x',num2str(1+(l-1)*nf)],num2str(xub,15)}),'^2'];
end
h_form{1} = [h_form{1},')'];
h_form = str2sym(h_form);
clhs = str2sym({'u1'});
crhs = cell(1,1);
crhs{1} = '(-2*(';
for l = 1:length(w)
    crhs{1} = [crhs{1},'+',num2str(w{l},15)];
    crhs{1} = [crhs{1},'*',replace(hl_form,{'s1','sub'},{['x',num2str(1+(l-1)*nf)],num2str(xub,15)})];
end
crhs{1} = [crhs{1},')*('];
for l = 1:length(w)
    crhs{1} = [crhs{1},'+',num2str(w{l},15)];
    crhs{1} = [crhs{1},'*',replace(dhldt0_form,{'s1','s2','th1'},{['x',num2str(1+(l-1)*nf)],['x',num2str(2+(l-1)*nf)],num2str(theta{l,1},15)})];
end
crhs{1} = [crhs{1},')+',num2str(r^2,15),'/',num2str(1+r^2,15),'*2*('];
for l = 1:length(w)
    crhs{1} = [crhs{1},'+',num2str(w{l},15)];
    crhs{1} = [crhs{1},'*',replace(hl_form,{'s1','sub'},{['x',num2str(1+(l-1)*nf)],num2str(xub,15)})];
    crhs{1} = [crhs{1},'*',replace(dhldt0_form,{'s1','s2','th1'},{['x',num2str(1+(l-1)*nf)],['x',num2str(2+(l-1)*nf)],num2str(theta{l,1},15)})];
end
crhs{1} = [crhs{1},'))/(-2*('];
for l = 1:length(w)
    crhs{1} = [crhs{1},'+',num2str(w{l},15)];
    crhs{1} = [crhs{1},'*',replace(hl_form,{'s1','sub'},{['x',num2str(1+(l-1)*nf)],num2str(xub,15)})];
end
crhs{1} = [crhs{1},')*('];
for l = 1:length(w)
    crhs{1} = [crhs{1},'+',num2str(w{l},15)];
    crhs{1} = [crhs{1},'*',replace(dhldt1_form,{'s1'},{['x',num2str(1+(l-1)*nf)]})];
end
crhs{1} = [crhs{1},')+',num2str(r^2,15),'/',num2str(1+r^2,15),'*2*('];
for l = 1:length(w)
    crhs{1} = [crhs{1},'+',num2str(w{l},15)];
    crhs{1} = [crhs{1},'*',replace(hl_form,{'s1','sub'},{['x',num2str(1+(l-1)*nf)],num2str(xub,15)})];
    crhs{1} = [crhs{1},'*',replace(dhldt1_form,{'s1'},{['x',num2str(1+(l-1)*nf)]})];
end
crhs{1} = [crhs{1},'))'];
crhs{1} = str2sym(crhs{1});
nci = 1;
end
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
deriv = 1;

[nci,phi,dphidxt,d2phidxt2,f,dfdxp,d2Hdxp2,d2fdxp2,h,dhdx] = hessian_compile(nci,clhs,crhs,lb,ub,phi_form,f_form,gh_form,h_form,np,ni,nt,nti,deriv,...
    {'theta1','theta2','Z','w';theta(:,1),theta(:,2),Z,w},[length(w),length(Z),length(w)],length(w));
nphi = size(phi_form,1);

idx_x0 = length(x0)+(1:ni);
av = [-1,1,-2];
ti = zeros(1,0);
tsc = [4.186312557738331,7.367536075801477];
tfc = 12;
x0c = [x0;0.736253127427679];
p0c = -0.165755283328793;
% av = [1,-1,-2];
% ti = zeros(1,0);
% tsc = [2.3197,12.0000-1e-5];
% tfc = 12;
% x0c = [x0;1.0000-1e-5];
% p0c = 0.0000;
% av = [1,-2,-1];
% ti = zeros(1,0);
% tsc = [5.3340,9.4772];
% tfc = 12;
% x0c = [x0;0.4817];
% p0c = -0.0903;
% av = [-1,1,-1];
% ti = zeros(1,0);
% tsc = [5.1300,10.1582];
% tfc = 12;
% x0c = [x0;1.0000-1e-5];
% p0c = 0.0000;
% av = [-1,1,-2];
% ti = zeros(1,0);
% tsc = [4.9775,12.0000-1e-5];
% tfc = 12;
% x0c = [x0;1.0000-1e-5];
% p0c = -0.1291;
% av = [-2,1,-1];
% ti = zeros(1,0);
% tsc = [1.7936,9.2595];
% tfc = 12;
% x0c = [x0;0.0000+1e-5];
% p0c = 0.1339;
% av = [-2,1,-2];
% ti = zeros(1,0);
% tsc = [2.3197,12.0000-1e-5];
% tfc = 12;
% x0c = [x0;0.0000+1e-5];
% p0c = 0.0000;
% av = [1,-1,-2];
% ti = zeros(1,0);
% tsc = [2.2512,12.0000-1e-5];
% tfc = 12;
% x0c = [x0;1.0000-1e-5];
% p0c = 0.0000;
% av = [1,-2,-1];
% ti = zeros(1,0);
% tsc = [5.2057,9.3303];
% tfc = 12;
% x0c = [x0;0.5378];
% p0c = -0.1029;
% av = [-1,1,-1];
% ti = zeros(1,0);
% tsc = [5.1155,10.3411];
% tfc = 12;
% x0c = [x0;1.0000-1e-5];
% p0c = 0.0000;
% av = [-1,1,-2];
% ti = zeros(1,0);
% tsc = [5.3782,12.0000-1e-5];
% tfc = 12;
% x0c = [x0;1.0000-1e-5];
% p0c = -0.0514;
% av = [-2,1,-1];
% ti = zeros(1,0);
% tsc = [1.5561,9.3240];
% tfc = 12;
% x0c = [x0;0.1044];
% p0c = 0.1153;
% av = [-2,1,-2];
% ti = zeros(1,0);
% tsc = [2.1886,12.0000-1e-5];
% tfc = 12;
% x0c = [x0;0.2359];
% p0c = -0.0240;
sc = [0.3,1,1,0.3,0.3];
hessian_check(nci,1:nphi,av,idx_x0,ti,tsc,tfc,x0c,p0c,sc,deriv);

tf = 12;
av = 3;%[1,-2,-1];
nt = 2;%3;
scv = [0.3,1,1,0.3,0.3];
deriv = 0;%5e-1/3;
n = 8;
m = 3000;
cheby = false;
suppl = 30;
[tau,p_sol,d_sol,deg,tel,fvalv,uv,N,p_mon,p_s,d,dus,flag,cphi,As,ys,y,Av,yv,u0,du,sc] = hessian_solve(tf,x0,nci,nphi,av,np,ni,nt,ti,scv,deriv,lb,ub,n,m,cheby,suppl);

profile viewer