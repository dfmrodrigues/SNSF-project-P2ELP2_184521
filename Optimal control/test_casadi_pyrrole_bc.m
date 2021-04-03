% An implementation of MPC using direct methods
% Adapted from Joel Andersson, 2016

clc
close all
clear
rng(42)

import casadi.*

format long

k1p = 0.053;
k2p = 0.128;
k3p = 0.028;
pp = [k1p;k2p;k3p];
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

t0 = 0;
tf = tfmax; % Time horizon
N = 50; % number of control intervals
cvode = true;
tv = linspace(t0,tf,N+1);

% Bounds on u
lbu = Fmin;
ubu = Fmax;

% Bounds on x
lbx = [-inf;-inf;-inf;-inf];
ubx = [inf;inf;inf;inf];

% Initial condition for x
x0 = [0.0;0.0;0.0;0.0];

% Declare model variables
t = SX.sym('t');
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');
x = [x1;x2;x3;x4];
u1 = SX.sym('u1');
u = u1;
k1 = SX.sym('k1');
k2 = SX.sym('k2');
k3 = SX.sym('k3');
p = [k1;k2;k3];
nx = size(x,1);
nu = size(u,1);
np = size(p,1);

% Model equations
nA = nA0+nAx{1}*x1+nAx{2}*x2+nAx{3}*x3+nAx{4}*x4;
nB = nB0+nBx{1}*x1+nBx{2}*x2+nBx{3}*x3+nBx{4}*x4;
nC = nC0+nCx{1}*x1+nCx{2}*x2+nCx{3}*x3+nCx{4}*x4;
nD = nD0+nDx{1}*x1+nDx{2}*x2+nDx{3}*x3+nDx{4}*x4;
V = V0+Vx{1}*x1+Vx{2}*x2+Vx{3}*x3+Vx{4}*x4;
xdot = [k1*nA*nB/V;k2*nB^2/V;k3*nB;u1/1e3];
k1m = 0.05;
k2m = 0.12;
k3m = 0.02;
pm = [k1m;k2m;k3m];
Spm = diag([0.01;0.025;0.005].^2);
Sigma = 0.05^2*eye(np);

% Objective term
qdot = 0*u1^2;
phi = -nC;
psi = [nB-cBmax*V;nD-cDmax*V];
npsi = size(psi,1);

% rng(42);
% H1 = rand(4,2);
% H2 = rand(2,2);
% H12 = [H1;H2];
% K = diag(rand(2,1));
% Sigma = rand(1)*eye(2);
% bg1 = (H1'*H1+Sigma/K)\([H1',zeros(2,2)]);
% Sigmag1 = Sigma/(H1'*H1+Sigma/K);
% bg2 = bg1+(H2'*H2+Sigma/Sigmag1)\([zeros(2,4),H2']-H2'*H2*bg1);
% Sigmag2 = Sigma/(H2'*H2+Sigma/Sigmag1);
% bg12 = (H12'*H12+Sigma/K)\(H12');
% Sigmag12 = Sigma/(H12'*H12+Sigma/K);
% bg2-bg12
% Sigmag2-Sigmag12

noise = true;
if(noise==true)
sigmaJ = 5e-3;
sigmaT = [1e-3;1e-3];
% invKJ = inv([0.018318157465348,-0.000010311687959,-0.000049160363310,-0.005757272907669,-0.747237299070200]'*...
%     [0.018318157465348,-0.000010311687959,-0.000049160363310,-0.005757272907669,-0.747237299070200]);
% invKT1 = inv([-0.010921879534378,-0.000019881352310,-0.000192675167824,-0.007468586184723,-1.224447430059715]'*...
%     [-0.010921879534378,-0.000019881352310,-0.000192675167824,-0.007468586184723,-1.224447430059715]);
% invKT2 = inv([-0.037643292272918,-0.000181183522937,-0.000122963008074,-0.043658112659854,-4.808230043366222]'*...
%     [-0.037643292272918,-0.000181183522937,-0.000122963008074,-0.043658112659854,-4.808230043366222]);
% invKJ = diag(([0,10,10,1,1e-3]/1e2).^2);
% invKT1 = diag(([0,10,10,1,1e-3]/1e2).^2);
% invKT2 = diag(([0,10,10,1,1e-3]/1e2).^2);
step = 1;%1e-2;
pp0 = pp;
pm0 = pm;
else
sigmaJ = 0;
sigmaT = [0;0];
% invKJ = [0,0,0,0,0];
% invKT1 = [0,0,0,0,0];
% invKT2 = [0,0,0,0,0];
step = 1e-6;
pp0 = pp;
pm0 = pp;
end
K_th = 1;
dt = 50;
niter = 16;
sc = [10,10,1,1e-3];
av = [-2,1,-1];
nt = length(av);
ni = sum(av>0);
ts0 = tf*(1:(nt-1))'/nt;
z00 = zeros(ni*nu,1);
p00 = zeros(ni*nu,1);
x0 = [x0;z00];
dJk = zeros(1,0);
dTk = zeros(npsi,0);
Uk = zeros(nt+2*ni*nu,0);
y = cell(3,niter+1);
S = cell(3,niter+1);
cm = cell(3,niter+1);
lambdam = cell(3,niter+1);
for l = 1:3
    y{l,1} = zeros(0,1);
    S{l,1} = zeros(0,0);
    cm{l,1} = 0;
    lambdam{l,1} = Inf;
end
u0 = [ts0;z00;p00];
eps0 = {y(:,1),S(:,1),cm(:,1),lambdam(:,1),u0,Uk,sc};
% eps_J0 = 0;
% eps_T0 = zeros(npsi,1);
% eps_dJdu0 = zeros(1,nt-1+2*ni*nu);
% eps_dTdu0 = zeros(npsi,nt-1+2*ni*nu);
% u0 = [ts0;z00;p00];
% eps0 = {eps_J0,eps_T0,eps_dJdu0,eps_dTdu0,u0,Uk,sc};
lam0 = {};
[tsp0,z0p0,p0p0,x0p0,lamp0,phip0,psip0,Fp0] =...
    direct_parsimonious(t0,tf,N,cvode,lbu,ubu,av,ts0,z00,p00,lbx,ubx,x0,t,x,u,p,pp0,xdot,qdot,phi,psi,eps0,lam0);
[tsm0,z0m0,p0m0,x0m0,lamm0,phim0,psim0,Fm0] =...
    direct_parsimonious(t0,tf,N,cvode,lbu,ubu,av,ts0,z00,p00,lbx,ubx,x0,t,x,u,p,pm0,xdot,qdot,phi,psi,eps0,lam0);
[J0,T0] = calc_parsimonious(tf,N,x0p0(1:nx),tsp0,z0p0,p0p0,pp,phip0,psip0,Fp0);
ts0 = tsm0;
z00 = z0m0;
p00 = p0m0;
Jpk = zeros(1,niter);
Tpk = zeros(npsi,niter);
% thJ = 0;
% thT = [0;0];
% sJpk = {};
% sTpk = {};
sigma{1} = sigmaJ(1);
sigma{2} = sigmaT(1);
sigma{3} = sigmaT(2);
for k = 1:niter
[Jp0,Tp0] = calc_parsimonious(tf,N,x0m0(1:nx),tsm0,z0m0,p0m0,pp,phim0,psim0,Fm0);
[Jm0,Tm0] = calc_parsimonious(tf,N,x0m0(1:nx),tsm0,z0m0,p0m0,pm,phim0,psim0,Fm0);
[Jppts1,Tppts1] = calc_parsimonious(tf,N,x0m0(1:nx),tsm0+[sc(1);0]*step,z0m0,p0m0,pp,phim0,psim0,Fm0);
[Jmpts1,Tmpts1] = calc_parsimonious(tf,N,x0m0(1:nx),tsm0+[sc(1);0]*step,z0m0,p0m0,pm,phim0,psim0,Fm0);
[Jppts2,Tppts2] = calc_parsimonious(tf,N,x0m0(1:nx),tsm0+[0;sc(2)]*step,z0m0,p0m0,pp,phim0,psim0,Fm0);
[Jmpts2,Tmpts2] = calc_parsimonious(tf,N,x0m0(1:nx),tsm0+[0;sc(2)]*step,z0m0,p0m0,pm,phim0,psim0,Fm0);
[Jppz01,Tppz01] = calc_parsimonious(tf,N,x0m0(1:nx),tsm0,z0m0+sc(3)*step,p0m0,pp,phim0,psim0,Fm0);
[Jmpz01,Tmpz01] = calc_parsimonious(tf,N,x0m0(1:nx),tsm0,z0m0+sc(3)*step,p0m0,pm,phim0,psim0,Fm0);
[Jppp01,Tppp01] = calc_parsimonious(tf,N,x0m0(1:nx),tsm0,z0m0,p0m0+sc(4)*step,pp,phim0,psim0,Fm0);
[Jmpp01,Tmpp01] = calc_parsimonious(tf,N,x0m0(1:nx),tsm0,z0m0,p0m0+sc(4)*step,pm,phim0,psim0,Fm0);
Jpk(:,k) = Jp0;
Tpk(:,k) = Tp0;
Jp0 = normrnd(Jp0,sigmaJ);
Tp0 = normrnd(Tp0,sigmaT);
Jppts1 = normrnd(Jppts1,sigmaJ);
Tppts1 = normrnd(Tppts1,sigmaT);
Jppts2 = normrnd(Jppts2,sigmaJ);
Tppts2 = normrnd(Tppts2,sigmaT);
Jppz01 = normrnd(Jppz01,sigmaJ);
Tppz01 = normrnd(Tppz01,sigmaT);
Jppp01 = normrnd(Jppp01,sigmaJ);
Tppp01 = normrnd(Tppp01,sigmaT);
dJk = [dJk(:,1:end),[Jp0,Jppts1,Jppts2,Jppz01,Jppp01]-[Jm0,Jmpts1,Jmpts2,Jmpz01,Jmpp01]];
dTk = [dTk(:,1:end),[Tp0,Tppts1,Tppts2,Tppz01,Tppp01]-[Tm0,Tmpts1,Tmpts2,Tmpz01,Tmpp01]];
Uk = [Uk(:,1:end),[1;tsm0-ts0;z0m0-z00;p0m0-p00]+[0,zeros(1,nt-1+2*ni*nu);zeros(nt-1+2*ni*nu,1),diag(sc*step)]];
dJk = dJk(:,max(1,end-dt+1):end);
dTk = dTk(:,max(1,end-dt+1):end);
Uk = Uk(:,max(1,end-dt+1):end);
y{1,k+1} = dJk(1,:)';
y{2,k+1} = dTk(1,:)';
y{3,k+1} = dTk(2,:)';
n = size(Uk,2);
In = eye(n);
sumsqs = 0;
sumsqv = sum(((Uk(:,2)-Uk)./[1;sc']).^2)';
sumsqm = cell2mat(arrayfun(@(i)sum(((Uk(:,i)-Uk)./[1;sc']).^2)',1:n,'UniformOutput',false));
[Lambda,C] = meshgrid(10.^(-1:0.02:1),10.^(-4:0.04:0));
K = cell(size(C,1),size(Lambda,2));
for i = 1:size(C,1)
    for j = 1:size(Lambda,2)
        K{i,j} = C(i,j)*exp(-sumsqm/2/Lambda(i,j)^2);
    end
end
Z = zeros(size(C,1),size(Lambda,2));
for l = 1:3
for i = 1:size(C,1)
    for j = 1:size(Lambda,2)
        Sij = K{i,j}+sigma{l}^2*In;
        Z(i,j) = (y{l,k+1}'/Sij*y{l,k+1}+log(det(Sij)))/n;
    end
end
cm{l,k+1} = min(C(Z==min(min(Z))));
lambdam{l,k+1} = min(Lambda(Z==min(min(Z))));
Ks = cm{l,k+1}*exp(-sumsqs/2/lambdam{l,k+1}^2);
Kv = cm{l,k+1}*exp(-sumsqv/2/lambdam{l,k+1}^2);
Km = cm{l,k+1}*exp(-sumsqm/2/lambdam{l,k+1}^2);
S{l,k+1} = Km+sigma{l}^2*In;
yb = Kv'/S{l,k+1}*y{l,k+1};
Syb = Ks-Kv'/S{l,k+1}*Kv;
yb0 = Km/S{l,k+1}*y{l,k+1};
Syb0 = Km-Km/S{l,k+1}*Km;
end
u0 = [ts0;z00;p00];
epsm0 = {y(:,k+1),S(:,k+1),cm(:,k+1),lambdam(:,k+1),u0,Uk,sc};
% thJ = (1-K_th)*thJ+K_th*dJk*Uk'/(Uk*Uk'+sigmaJ^2*invKJ);
% thT = [(1-K_th)*thT(1)+K_th*dTk(1,:)*Uk'/(Uk*Uk'+sigmaT(1)^2*invKT1);...
%     (1-K_th)*thT(2)+K_th*dTk(2,:)*Uk'/(Uk*Uk'+sigmaT(2)^2*invKT2)];
% sJpk = [sJpk,sigmaJ^2*eye(5)/(Uk*Uk'+sigmaJ^2*invKJ)];
% sTpk = [sTpk,{sigmaT(1)^2*eye(5)/(Uk*Uk'+sigmaT(1)^2*invKT1);...
%     sigmaT(2)^2*eye(5)/(Uk*Uk'+sigmaT(2)^2*invKT2)}];
% eps_Jm0 = thJ(:,1);
% eps_Tm0 = thT(:,1);
% eps_dJdum0 = thJ(:,2:end);
% eps_dTdum0 = thT(:,2:end);
% u0 = [ts0;z00;p00];
% epsm0 = {eps_Jm0,eps_Tm0,eps_dJdum0,eps_dTdum0,u0,Uk,sc};
[tsm0,z0m0,p0m0,x0m0,lamm0,phim0,psim0,Fm0] =...
    direct_parsimonious(t0,tf,N,cvode,lbu,ubu,av,tsm0,z0m0,p0m0,lbx,ubx,x0,t,x,u,p,pm0,xdot,qdot,phi,psi,epsm0,lam0);
end
figure,hold on;
yyaxis left;
plot(1:niter,Jpk','-','LineWidth',1.5);
plot(1:niter,J0*ones(1,niter)','--','LineWidth',1.5);
xlabel('$k$','Interpreter','LaTeX','FontSize',22);
ylabel('$\hat{\phi}^{p}(${\boldmath$\tau$}$_{k})$','Interpreter','LaTeX','FontSize',22);
set(gca,'FontSize',20,'XLim',[1,16],'XTick',1:3:16);
yyaxis right;
plot(1:niter,Tpk','-','LineWidth',1.5);
plot(1:niter,zeros(1,niter)','--','LineWidth',1.5);
ylabel('$\hat{\psi_{j}}^{p}(${\boldmath$\tau$}$_{k})$','Interpreter','LaTeX','FontSize',22);
set(gca,'Position',[0.18,0.135,0.64,0.78]);
set(gca,'FontSize',20,'XLim',[1,16],'XTick',1:3:16);
figure,plot(cell2mat(cm'))
figure,plot(cell2mat(lambdam'))
% figure,plot(cellfun(@(e)e(1,1),sJpk)')
% figure,plot(cellfun(@(e)e(1,1),sTpk)')
keyboard
pm = pp;

% Initial solution using direct parsimonious input parameterization
pmp = pm;
Spmp = Spm;
av = [-2,1,-1];
nt = length(av);
ni = sum(av>0);
ts0 = tf*(1:(nt-1))'/nt;
z00 = zeros(ni*nu,1);
p00 = zeros(ni*nu,1);
x0 = [x0;z00];
eps_J0 = 0;
eps_T0 = zeros(npsi,1);
eps_dJdu0 = zeros(1,nt-1+2*ni*nu);
eps_dTdu0 = zeros(npsi,nt-1+2*ni*nu);
u0 = [ts0;z00;p00];
eps0 = {eps_J0,eps_T0,eps_dJdu0,eps_dTdu0,u0};
lam0 = {};
[ts_optp,z0_optp,p0_optp,x0_optp,lam_optp,phi_optp,psi_optp,F_optp] =...
    direct_parsimonious(t0,tf,N,cvode,lbu,ubu,av,ts0,z00,p00,lbx,ubx,x0,t,x,u,p,pmp,xdot,qdot,phi,psi,eps0,lam0);
[J_optp,T_optp,x_optp] = calc_parsimonious(tf,N,x0_optp(1:nx),ts_optp,z0_optp,p0_optp,pmp,phi_optp,psi_optp,F_optp);

% Simulate MPC
J_simp = zeros(1,N+1);
T_simp = zeros(npsi,N+1);
ts_simp = zeros(nt-1,N+1);
z0_simp = zeros(ni*nu,N+1);
p0_simp = zeros(ni*nu,N+1);
x_simp = zeros(size(x_optp));
ts0k = ts_optp;
z00k = z0_optp;
p00k = p0_optp;
x0k = x_optp(:,1);
lam0k = lam_optp;
for k = 1:(N+1)
    x_simp(:,k) = x0k(:,1);
    t0k = tv(k);
    [ts_optk,z0_optk,p0_optk,x0_optk,lam_optk,phi_optk,psi_optk,F_optk] =...
        direct_parsimonious(t0k,tf,N-k+1,cvode,lbu,ubu,av,ts0k,z00k,p00k,lbx,ubx,x0k,t,x,u,p,pmp,xdot,qdot,phi,psi,eps0,lam0k);
    [J_optk,T_optk,x_optk] = calc_parsimonious(tf,N-k+1,x0_optk(1:nx),ts_optk,x0_optk((nx+1):end),p0_optk,pmp,phi_optk,psi_optk,F_optk);
    disp(sum(sum(abs(ts_optk-ts_optp))))
    disp(sum(sum(abs(z0_optk-z0_optp))))
    disp(sum(sum(abs(p0_optk-p0_optp))))
    disp(sum(sum(abs(x_optk(:,1:end)-x_optp(:,k:end)))))
    J_simp(:,k) = J_optk;
    T_simp(:,k) = T_optk;
    ts_simp(:,k) = ts_optk;
    z0_simp(:,k) = z0_optk;
    p0_simp(:,k) = p0_optk;
    if(k==N+1)
        break;
    end
    ts0k = ts_optk;
    z00k = z0_optk;
    p00k = p0_optk;
    fk = F_optk{1}('x0',x_optk(:,1),'p',pp,'ts',ts_optk,'z0',z0_optk,'p0',p0_optk);
    x0k = full(fk.xf);
    lam0k = lam_optk;
end
disp(ts_simp-ts_optp)
disp(z0_simp-z0_optp)
disp(p0_simp-p0_optp)
disp(x_simp-x_optp)

% Plot the solution
for c = 'os'
if(c=='o')
    ts_opt = ts_optp;
    z0_opt = z0_optp;
    p0_opt = p0_optp;
    x_opt = x_optp;
elseif(c=='s')
    ts_opt = ts_simp(:,end);
    z0_opt = z0_simp(:,end);
    p0_opt = p0_simp(:,end);
    x_opt = x_simp;
end
x1_opt = x_opt(1,:);
x2_opt = x_opt(2,:);
x3_opt = x_opt(3,:);
x4_opt = x_opt(4,:);
figure;
hold on
plot(tv,x1_opt,tv,x2_opt,tv,x3_opt,tv,x4_opt)
teb = [t0,ts_opt(1:(nt-1))',tf];
prev = [];
for i = 1:nt
    if(av(i)>0)
        hold on;ax = gca;ax.ColorOrderIndex = 5;
        j = av(i);
        plot([teb(i),teb(i)],[prev,z0_opt(j)]);
        hold on;ax = gca;ax.ColorOrderIndex = 5;
        prev = z0_opt(j)+p0_opt(j)*(teb(i+1)-teb(i));
        plot([teb(i),teb(i+1)],[z0_opt(j),prev]);
    elseif(av(i)==-1)
        hold on;ax = gca;ax.ColorOrderIndex = 5;
        plot([teb(i),teb(i)],[prev,lbu(1)]);
        hold on;ax = gca;ax.ColorOrderIndex = 5;
        prev = lbu(1);
        plot([teb(i),teb(i+1)],[lbu(1),prev]);
    elseif(av(i)==-2)
        hold on;ax = gca;ax.ColorOrderIndex = 5;
        plot([teb(i),teb(i)],[prev,ubu(1)]);
        hold on;ax = gca;ax.ColorOrderIndex = 5;
        prev = ubu(1);
        plot([teb(i),teb(i+1)],[ubu(1),prev]);
    end
end
xlabel('t')
legend('x1','x2','x3','x4','u')
end
x1_simp = x_simp(1,:);
x2_simp = x_simp(2,:);
x3_simp = x_simp(3,:);
x4_simp = x_simp(4,:);
nA_simp = nA0+nAx{1}*x1_simp+nAx{2}*x2_simp+nAx{3}*x3_simp+nAx{4}*x4_simp;
nB_simp = nB0+nBx{1}*x1_simp+nBx{2}*x2_simp+nBx{3}*x3_simp+nBx{4}*x4_simp;
nC_simp = nC0+nCx{1}*x1_simp+nCx{2}*x2_simp+nCx{3}*x3_simp+nCx{4}*x4_simp;
nD_simp = nD0+nDx{1}*x1_simp+nDx{2}*x2_simp+nDx{3}*x3_simp+nDx{4}*x4_simp;
V_simp = V0+Vx{1}*x1_simp+Vx{2}*x2_simp+Vx{3}*x3_simp+Vx{4}*x4_simp;
xrdot_p_simp = [nA_simp.*nB_simp./V_simp;nB_simp.^2./V_simp;nB_simp];
y_simp = x_simp(1:3,end);
H_simp = diag(sum(xrdot_p_simp(:,1:end-1)+xrdot_p_simp(:,2:end),2)/2*(tf-t0)/N);
pmp = pmp+(H_simp'*H_simp+Sigma/Spmp)\(H_simp'*(y_simp-H_simp*pmp));
Spmp = Sigma/(H_simp'*H_simp+Sigma/Spmp);
keyboard

% Initial solution using direct single shooting
pms = pm;
Spms = Spm;
u0 = zeros(nu,N);
x0 = x0(1:nx);
lam0 = {};
[u_opts,x0_opts,lam_opts,phi_opts,psi_opts,F_opts] =...
    direct_single_shooting(t0,tf,N,cvode,lbu,ubu,u0,lbx,ubx,x0,t,x,u,p,pms,xdot,qdot,phi,psi,lam0);
[J_opts,T_opts,x_opts] = calc_single_shooting(tf,N,x0_opts,u_opts,pms,phi_opts,psi_opts,F_opts);

% Simulate MPC
J_sims = zeros(1,N+1);
T_sims = zeros(npsi,N+1);
u_sims = zeros(size(u_opts));
x_sims = zeros(size(x_opts));
u0k = u_opts;
x0k = x_opts(:,1);
lam0k = lam_opts;
for k = 1:(N+1)
    x_sims(:,k) = x0k(:,1);
    t0k = tv(k);
    [u_optk,x0_optk,lam_optk,phi_optk,psi_optk,F_optk] =...
        direct_single_shooting(t0k,tf,N-k+1,cvode,lbu,ubu,u0k,lbx,ubx,x0k,t,x,u,p,pms,xdot,qdot,phi,psi,lam0k);
    [J_optk,T_optk,x_optk] = calc_single_shooting(tf,N-k+1,x0_optk,u_optk,pms,phi_optk,psi_optk,F_optk);
    disp(sum(sum(abs(u_optk(:,1:end-1)-u_opts(:,k:end-1)))))
    disp(sum(sum(abs(x_optk(:,1:end)-x_opts(:,k:end)))))
    J_sims(:,k) = J_optk;
    T_sims(:,k) = T_optk;
    u_sims(:,k) = u_optk(:,1);
    if(k==N+1)
        break;
    end
    u0k = u_optk(:,2:end);
    fk = F_optk{1}('x0',x_optk(:,1),'p',pp,'u',u_optk(:,1));
    x0k = full(fk.xf);
    lam_u_optk = lam_optk{1};
    lam_x_optk = lam_optk{2};
    lam_dx_optk = lam_optk{3};
    lam_T_optk = lam_optk{4};
    lam0k = {lam_u_optk(:,2:end),lam_x_optk(:,2:end),lam_dx_optk(:,2:end),lam_T_optk};
end
disp(u_sims-u_opts)
disp(x_sims-x_opts)

% Plot the solution
for c = 'os'
if(c=='o')
    u_opt = u_opts;
    x_opt = x_opts;
elseif(c=='s')
    u_opt = u_sims;
    x_opt = x_sims;
end
x1_opt = x_opt(1,:);
x2_opt = x_opt(2,:);
x3_opt = x_opt(3,:);
x4_opt = x_opt(4,:);
u1_opt = u_opt(1,:);
figure;
hold on
plot(tv,x1_opt,tv,x2_opt,tv,x3_opt,tv,x4_opt)
stairs(tv,u1_opt)
xlabel('t')
legend('x1','x2','x3','x4','u')
end
x1_sims = x_sims(1,:);
x2_sims = x_sims(2,:);
x3_sims = x_sims(3,:);
x4_sims = x_sims(4,:);
nA_sims = nA0+nAx{1}*x1_sims+nAx{2}*x2_sims+nAx{3}*x3_sims+nAx{4}*x4_sims;
nB_sims = nB0+nBx{1}*x1_sims+nBx{2}*x2_sims+nBx{3}*x3_sims+nBx{4}*x4_sims;
nC_sims = nC0+nCx{1}*x1_sims+nCx{2}*x2_sims+nCx{3}*x3_sims+nCx{4}*x4_sims;
nD_sims = nD0+nDx{1}*x1_sims+nDx{2}*x2_sims+nDx{3}*x3_sims+nDx{4}*x4_sims;
V_sims = V0+Vx{1}*x1_sims+Vx{2}*x2_sims+Vx{3}*x3_sims+Vx{4}*x4_sims;
xrdot_p_sims = [nA_sims.*nB_sims./V_sims;nB_sims.^2./V_sims;nB_sims];
y_sims = x_sims(1:3,end);
H_sims = diag(sum(xrdot_p_sims(:,1:end-1)+xrdot_p_sims(:,2:end),2)/2*(tf-t0)/N);
pms = pms+(H_sims'*H_sims+Sigma/Spms)\(H_sims'*(y_sims-H_sims*pms));
Spms = Sigma/(H_sims'*H_sims+Sigma/Spms);

% Initial solution using direct multiple shooting
pmm = pm;
Spmm = Spm;
u0 = zeros(nu,N);
x0 = x0(1:nx);
lam0 = {};
[u_optm,x0_optm,lam_optm,phi_optm,psi_optm,F_optm] =...
    direct_multiple_shooting(t0,tf,N,cvode,lbu,ubu,u0,lbx,ubx,x0,t,x,u,p,pmm,xdot,qdot,phi,psi,lam0);
[J_optm,T_optm,x_optm] = calc_multiple_shooting(tf,N,x0_optm,u_optm,pmm,phi_optm,psi_optm,F_optm);

% Simulate MPC
J_simm = zeros(1,N+1);
T_simm = zeros(npsi,N+1);
u_simm = zeros(size(u_optm));
x_simm = zeros(size(x_optm));
u0k = u_optm;
x0k = x_optm(:,1);
lam0k = lam_optm;
for k = 1:(N+1)
    x_simm(:,k) = x0k(:,1);
    t0k = tv(k);
    [u_optk,x0_optk,lam_optk,phi_optk,psi_optk,F_optk] =...
        direct_multiple_shooting(t0k,tf,N-k+1,cvode,lbu,ubu,u0k,lbx,ubx,x0k,t,x,u,p,pmm,xdot,qdot,phi,psi,lam0k);
    [J_optk,T_optk,x_optk] = calc_multiple_shooting(tf,N-k+1,x0_optk,u_optk,pmm,phi_optk,psi_optk,F_optk);
    disp(sum(sum(abs(u_optk(:,1:end-1)-u_optm(:,k:end-1)))))
    disp(sum(sum(abs(x_optk(:,1:end)-x_optm(:,k:end)))))
    J_simm(:,k) = J_optk;
    T_simm(:,k) = T_optk;
    u_simm(:,k) = u_optk(:,1);
    if(k==N+1)
        break;
    end
    u0k = u_optk(:,2:end);
    fk = F_optk{1}('x0',x_optk(:,1),'p',pp,'u',u_optk(:,1));
    x0k = full(fk.xf);
    lam_u_optk = lam_optk{1};
    lam_x_optk = lam_optk{2};
    lam_dx_optk = lam_optk{3};
    lam_T_optk = lam_optk{4};
    lam0k = {lam_u_optk(:,2:end),lam_x_optk(:,2:end),lam_dx_optk(:,2:end),lam_T_optk};
end
disp(u_simm-u_optm)
disp(x_simm-x_optm)

% Plot the solution
for c = 'os'
if(c=='o')
    u_opt = u_optm;
    x_opt = x_optm;
elseif(c=='s')
    u_opt = u_simm;
    x_opt = x_simm;
end
x1_opt = x_opt(1,:);
x2_opt = x_opt(2,:);
x3_opt = x_opt(3,:);
x4_opt = x_opt(4,:);
u1_opt = u_opt(1,:);
figure;
hold on
plot(tv,x1_opt,tv,x2_opt,tv,x3_opt,tv,x4_opt)
stairs(tv,u1_opt)
xlabel('t')
legend('x1','x2','x3','x4','u')
end
x1_simm = x_simm(1,:);
x2_simm = x_simm(2,:);
x3_simm = x_simm(3,:);
x4_simm = x_simm(4,:);
nA_simm = nA0+nAx{1}*x1_simm+nAx{2}*x2_simm+nAx{3}*x3_simm+nAx{4}*x4_simm;
nB_simm = nB0+nBx{1}*x1_simm+nBx{2}*x2_simm+nBx{3}*x3_simm+nBx{4}*x4_simm;
nC_simm = nC0+nCx{1}*x1_simm+nCx{2}*x2_simm+nCx{3}*x3_simm+nCx{4}*x4_simm;
nD_simm = nD0+nDx{1}*x1_simm+nDx{2}*x2_simm+nDx{3}*x3_simm+nDx{4}*x4_simm;
V_simm = V0+Vx{1}*x1_simm+Vx{2}*x2_simm+Vx{3}*x3_simm+Vx{4}*x4_simm;
xrdot_p_simm = [nA_simm.*nB_simm./V_simm;nB_simm.^2./V_simm;nB_simm];
y_simm = x_simm(1:3,end);
H_simm = diag(sum(xrdot_p_simm(:,1:end-1)+xrdot_p_simm(:,2:end),2)/2*(tf-t0)/N);
pmm = pmm+(H_simm'*H_simm+Sigma/Spmm)\(H_simm'*(y_simm-H_simm*pmm));
Spmm = Sigma/(H_simm'*H_simm+Sigma/Spmm);
