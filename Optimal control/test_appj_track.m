clc
close all
clear
profile on

import casadi.*

tfmax = 120;
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
k1p = m2*eta/Rho/Cp/pi/r^2/d;%*0.8;
k2p = 2*pi*r*d*K*m1/Rho/Cp/pi/r^2/d*(1/2+(Tmax-Tb)/(Tb-Tinf))*log(1+(Tb-Tinf)/(Tmax-Tb));%0;
k3p = 0;%2*pi*r*d*K*m1/Rho/Cp/pi/r^2/d*0.8;
pp = [k1p;k2p;k3p];

t0 = 0;
tf = tfmax; % Time horizon
N = 50; % number of control intervals
cvode = true;

% Bounds on u
lbu = Pmin;
ubu = Pmax;

% Bounds on x
lbx = [-inf;-inf];
ubx = [inf;Tmax];

% Initial condition for x
x0 = [0.0;T0];

% Declare model variables
t = SX.sym('t');
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];
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
CEM = x1;
T = x2;
P = u1;
xdot = [0.5^(43.0+273.15-T)/60.0;k1*P-k2*((T-Tinf)-log(1+exp(10*(T-Tb)))/10)./(log(T-Tinf)-log(log(1+exp(10*(T-Tb)))/10))-k3*(T-(Tinf+Tb)/2)];
% xdot = [0.5^(43.0+273.15-T)/60.0;k1*P-k3*(T-(Tinf+Tb)/2)];
k1m = m2*eta/Rho/Cp/pi/r^2/d*0.8;
k2m = 0;
k3m = 2*pi*r*d*K*m1/Rho/Cp/pi/r^2/d*0.8;
pm = [k1m;k2m;k3m];
Spm = diag([0.01;0.025;0.005].^2);

% Objective term
qdot = 0*u1^2;
phi = t;
psi = [CEMsp-CEM;T-T0];
npsi = size(psi,1);

% Initial solution using direct single shooting
pms = pm;
Spms = Spm;
tf0 = tf;
u0 = (lbu+ubu)/2*ones(1,N);
x0 = x0(1:nx);
lam0 = {};
[tf_opts,u_opts,x0_opts,lam_opts,phi_opts,psi_opts,F_opts] =...
    direct_single_shooting(t0,tf,N,cvode,lbu,ubu,tf0,u0,lbx,ubx,x0,t,x,u,p,pms,xdot,qdot,phi,psi,lam0);
tv_opts = linspace(t0,tf_opts,N+1);
[J_opts,T_opts,x_opts] = calc_single_shooting(N,x0_opts,tf_opts,u_opts,pms,phi_opts,psi_opts,F_opts);

% Simulate MPC
J_sims = zeros(1,N+1);
T_sims = zeros(npsi,N+1);
tf_sims = zeros(size(tf_opts));
u_sims = zeros(size(u_opts));
tv_sims = zeros(1,N+1);
x_sims = zeros(size(x_opts));
tf0k = tf_opts;
u0k = u_opts;
x0k = x_opts(:,1);
t0k = t0;
lam0k = lam_opts;
for k = 1:(N+1)
    x_sims(:,k) = x0k(:,1);
    tv_sims(k) = t0k;
    [tf_optk,u_optk,x0_optk,lam_optk,phi_optk,psi_optk,F_optk] =...
        direct_single_shooting(t0k,tf,N-k+1,cvode,lbu,ubu,tf0k,u0k,lbx,ubx,x0k,t,x,u,p,pms,xdot,qdot,phi,psi,lam0k);
    [J_optk,T_optk,x_optk] = calc_single_shooting(N-k+1,x0_optk,tf_optk,u_optk,pms,phi_optk,psi_optk,F_optk);
    disp(sum(sum(abs(tf_optk-tf_opts))))
    disp(sum(sum(abs(u_optk(:,1:end-1)-u_opts(:,k:end-1)))))
    disp(sum(sum(abs(x_optk(:,1:end)-x_opts(:,k:end)))))
    J_sims(:,k) = J_optk;
    T_sims(:,k) = T_optk;
    tf_sims(:,k) = tf_optk;
    u_sims(:,k) = u_optk(:,1);
    if(k==N+1)
        break;
    end
    tf0k = tf_optk;
    u0k = u_optk(:,2:end);
    fk = F_optk{1}('x0',x_optk(:,1),'p',pp,'tf',tf_optk,'u',u_optk(:,1));
    x0k = full(fk.xf);
    t0k = t0k+(tf_optk-t0k)/(N-k+1);
    lam_tf_optk = lam_optk{1};
    lam_u_optk = lam_optk{2};
    lam_dt_optk = lam_optk{3};
    lam_x_optk = lam_optk{4};
    lam_dx_optk = lam_optk{5};
    lam_T_optk = lam_optk{6};
    lam0k = {lam_tf_optk,lam_u_optk(:,2:end),lam_dt_optk,lam_x_optk(:,2:end),lam_dx_optk(:,2:end),lam_T_optk};
end
disp(tf_sims-tf_opts)
disp(u_sims-u_opts)
disp(x_sims-x_opts)

% Plot the solution
tP = u_sims;
tout = tv_sims;
CEM = x_sims(1,:);
T = x_sims(2,:);
Tr0 = @(t)Tmax*ones(size(t));
figure;
ax{1} = subplot(1,3,1);
ax{2} = subplot(1,3,2);
ax{3} = subplot(1,3,3);
hold(ax{1},'on'),plot(ax{1},tout,T-T0,'b','LineWidth',2);set(ax{1},'Children',flipud(get(ax{1},'Children')));
hold(ax{1},'on'),plot(ax{1},tout,Tr0(tout)-T0,'r');
hold(ax{1},'on'),plot(ax{1},tout,zeros(size(tout)),'r--');
title(ax{1},'');ylabel(ax{1},'$T(t)-T_{0}$ [K]','Interpreter','LaTeX','FontSize',20);
p = xlabel(ax{1},'');set(p,'String','$t$ [s]','Interpreter','LaTeX','FontSize',20);
set(ax{1},'FontSize',18,'XTick',0:30:90,'XTickLabel',0:30:90);
hold(ax{2},'on'),plot(ax{2},tout,CEM,'b','LineWidth',2);set(ax{2},'Children',flipud(get(ax{2},'Children')));
hold(ax{2},'on'),plot(ax{2},tout,CEMsp*ones(size(tout)),'r--');
title(ax{2},'');ylabel(ax{2},'$CEM(t)$ [min]','Interpreter','LaTeX','FontSize',20);
p = xlabel(ax{2},'');set(p,'String','$t$ [s]','Interpreter','LaTeX','FontSize',20);
set(ax{2},'FontSize',18,'XTick',0:30:90,'XTickLabel',0:30:90);
hold(ax{3},'on'),stairs(ax{3},tout,tP(1,:),'b','LineWidth',2);set(ax{3},'Children',flipud(get(ax{3},'Children')));
hold(ax{3},'on'),plot(ax{3},tout,Pmin*ones(size(tout)),'r');
hold(ax{3},'on'),plot(ax{3},tout,Pmax*ones(size(tout)),'r');
title(ax{3},'');ylabel(ax{3},'$\tilde{P}(t)$ [W]','Interpreter','LaTeX','FontSize',20);
p = xlabel(ax{3},'');set(p,'String','$t$ [s]','Interpreter','LaTeX','FontSize',20);
set(ax{3},'FontSize',18,'XTick',0:30:90,'XTickLabel',0:30:90);
set(gcf,'Units','normalized'),set(gcf,'OuterPosition',[1,1,2,1].*get(gcf,'OuterPosition'));
keyboard

clc
close all
clear

% Define constants for the model
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
T0 = 310.15;%312.15; % K
Pmin = 1.0; % W
Pmax = 5.0; % W
Tmax = 316.15; % K
CEMsp = 1.5; % min
k1 = m2*eta/Rho/Cp/pi/r^2/d*0.8;
k2 = 2*pi*r*d*K*m1/Rho/Cp/pi/r^2/d*0.8;

% Write dynamic model of APPJ
glhs = str2sym({'CEM','T','P'}); % Common name of variables
grhs = str2sym({'x1','x2','u1'}); % Name of variables for symbolic manipulation
phi_form = subs(str2sym({'t/100';'T-T0';'CEMsp-CEM'}),glhs,grhs); % Terminal cost and constraints
f_form = subs(str2sym({'0.5^(43.0+273.15-T)/60.0';'k1*P-k2*(T-(Tinf+Tb)/2)'}),glhs,grhs); % Differential equations
gh_form = str2sym({});
h_form = subs(str2sym({'T-Tmax'}),glhs,grhs); % State path constraints
x0 = str2sym({'0.0';'T0'}); % Initial conditions
lb = str2sym({'Pmin'}); % Lower input bound
ub = str2sym({'Pmax'}); % Upper input bound
clhs = str2sym({'u1'}); % Left-hand side of input assignment in path constraints
crhs = {subs(str2sym({'k2*(T-(Tinf+Tb)/2)/k1'}),glhs,grhs)}; % Right-hand side of input assignment in path constraints
nci = 1; % One case of input assignment
np = 0; % No free parameters
ni = 0; % No initial condition for sensitivity-seeking arcs
nt = 2; % Two switching times, of which one final time
nti = 0; % No intermediate times considered in cost and constraints
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

% Compile dynamic model for efficient numerical solution of optimal control problem (OCP)
[nci,phi,dphidxt,d2phidxt2,f,dfdxp,d2Hdxp2,d2fdxp2,h,dhdx] = hessian_compile(nci,clhs,crhs,lb,ub,phi_form,f_form,gh_form,h_form,np,ni,nt,nti,deriv);
nphi = length(phi_form);

% Compile auxiliary files for numerical integrators
eval('mex -outdir ../ode45_rp -output rp45_mex ../ode45_rp/rp45_mex.c ../ode45_rp/rp45.c ../ode45_rp/ntrp45.c');
eval('mex -outdir ../ode45_rp -output ntrp45_mex ../ode45_rp/ntrp45_mex.c ../ode45_rp/ntrp45.c');
eval('mex -outdir ../ode45_rp -output binterp_mex ../ode45_rp/binterp_mex.c ../ode45_rp/binterp.c');
eval('mex -outdir ../ode45_rp -output rp3h_mex ../ode45_rp/rp3h_mex.c ../ode45_rp/rp3h.c ../ode45_rp/ntrp3h.c');
eval('mex -outdir ../ode45_rp -output ntrp3h_mex ../ode45_rp/ntrp3h_mex.c ../ode45_rp/ntrp3h.c');
addpath('../ode45_rp');

% Check whether the gradient computation works well (it does)
idx_x0 = length(x0)+(1:ni);
av = [-2,-1];
ti = zeros(1,0);
tsc = 92.8404967453808;%91.178807262324270;
tfc = 109.4495929298785;%99.156543394002369;
x0c = x0;
p0c = [];
sc = [100,100];
hessian_check(nci,1:nphi,av,idx_x0,ti,tsc,tfc,x0c,p0c,sc,deriv);

% Determine the optimal solution for the model, which is the initial point for modifier adaptation
tf = 120; % Upper bound for the final time
av = [-2,-1]; % Arc sequence (input at upper bound, input at lower bound)
nt = 2; % Two switching times, of which one final time
scv = [100,100]; % Scaling of decision variables
deriv = 1e-2; % Scaling of slope in sensitivity-seeking arcs (not used)
% Commented code below would be used for global solution to the OCP
% n = 6;
% m = 1000;
% cheby = false;
% suppl = 0;
% xunits = 's';
% yunits = 'W';
% [tau,p_sol,d_sol,deg,tel,fvalv,uv,N,p_mon,p_s,d,dus,flag,cphi,As,ys,y,Av,yv,u0,du,~] = hessian_solve(tf,x0,nci,nphi,av,np,ni,nt,ti,scv,deriv,lb,ub,n,m,cheby,suppl,xunits,yunits);
% u = uv{1};
% Initialization of random number seed to make the random number sequence equal to the case of global solution to the OCP
rng(42,'twister');
for i = 1:1200
    rand(1,2);
end
rand(2,1);
rand(2,1);
lam0 = zeros(nphi,nt+length(x0)+1); % Initial values of modifiers equal zero
pi00 = zeros(nt+length(x0),1); % Initial values of decision variables not used
u = hessian_fmincon(nci,nphi,av,np,idx_x0,ti,tsc,tfc,x0,[],sc,lam0,pi00,lb,ub,0.001); % Local solution to the OCP for the unmodified model
tsm0 = u(1); % Optimal switching time to arc with input at lower bound for the unmodified model
tfm0 = u(2); % Optimal final time for the unmodified model

profile viewer

% Implement modifier adaptation
rng(42);
ts0 = tsm0; % Initial value of switching time to arc with input at lower bound
tf0 = tfm0; % Initial value of final time
step = 1e-2; % Step away from nominal point to compute gradients via finite differences is step*sc
npsi = 2; % Number of terminal constraints
dt = 50; % Number of previous points for computation of modifiers
niter = 5; % Number of iterations of modifier adaptation
dJk = zeros(1,0);
dTk = zeros(npsi,0);
Uk = zeros(nt+1,0);
Jpk = zeros(1,niter);
Tpk = zeros(npsi,niter);
for k = 1:niter
addpath('./hessian');
% Compute cost and constraints of the plant for the nominal point
[Tfp0,CEMfp0] = appj_track(tsm0,tfm0);
Jp0 = tfm0;
Tp0 = [Tfp0;CEMfp0];
% Compute cost and constraints of the modified model for the nominal point
Tfm0 = hessian_calc(nci,2,av,ti,0,tsm0,tfm0,x0,[],lam0,pi00);
CEMfm0 = hessian_calc(nci,3,av,ti,0,tsm0,tfm0,x0,[],lam0,pi00);
Jm0 = tfm0;
Tm0 = cell2mat([Tfm0;CEMfm0]);
% Compute cost and constraints of the plant for the point with perturbed switching time
[Tfppts1,CEMfppts1] = appj_track(tsm0+sc(1)*step,tfm0);
Jppts1 = tfm0;
Tppts1 = [Tfppts1;CEMfppts1];
% Compute cost and constraints of the modified model for the point with perturbed switching time
Tfmpts1 = hessian_calc(nci,2,av,ti,0,tsm0+sc(1)*step,tfm0,x0,[],lam0,pi00);
CEMfmpts1 = hessian_calc(nci,3,av,ti,0,tsm0+sc(1)*step,tfm0,x0,[],lam0,pi00);
Jmpts1 = tfm0;
Tmpts1 = cell2mat([Tfmpts1;CEMfmpts1]);
% Compute cost and constraints of the plant for the point with perturbed final time
[Tfpptf,CEMfpptf] = appj_track(tsm0,tfm0+sc(2)*step);
Jpptf = tfm0+sc(2)*step;
Tpptf = [Tfpptf;CEMfpptf];
% Compute cost and constraints of the modified model for the point with perturbed final time
Tfmptf = hessian_calc(nci,2,av,ti,0,tsm0,tfm0+sc(2)*step,x0,[],lam0,pi00);
CEMfmptf = hessian_calc(nci,3,av,ti,0,tsm0,tfm0+sc(2)*step,x0,[],lam0,pi00);
Jmptf = tfm0+sc(2)*step;
Tmptf = cell2mat([Tfmptf;CEMfmptf]);
% Save values of plant cost and constraints for plots
Jpk(:,k) = Jp0;
Tpk(:,k) = Tp0;
% Append plant-model mismatch in the cost and constraints to the relevant matrices
dJk = [dJk(:,1:end),[Jp0,Jppts1,Jpptf]-[Jm0,Jmpts1,Jmptf]];
dTk = [dTk(:,1:end),[Tp0,Tppts1,Tpptf]-[Tm0,Tmpts1,Tmptf]];
% Append distance to the initial values of the decision variables to the relevant matrix
Uk = [Uk(:,1:end),[1;tsm0-ts0;tfm0-tf0]+[0,zeros(1,nt);zeros(nt,1),diag(sc*step)]];
% Exclude points that were visited much before (not used in this case)
dJk = dJk(:,max(1,end-dt+1):end);
dTk = dTk(:,max(1,end-dt+1):end);
Uk = Uk(:,max(1,end-dt+1):end);
lam = [[dJk*Uk'/(Uk*Uk');dTk(1,:)*Uk'/(Uk*Uk');dTk(2,:)*Uk'/(Uk*Uk')],zeros(nphi,length(x0))]; % Compute updated value of modifiers
pi0 = [ts0;tf0;x0]; % Compute initial values of decision variables
% Commented code below would be used for global solution to the OCP
% [tau,p_sol,d_sol,deg,tel,fvalv,uv,N,p_mon,p_s,d,dus,flag,cphi,As,ys,y,Av,yv,u0,du,~] = hessian_solve(tf,x0,nci,nphi,av,np,ni,nt,ti,scv,deriv,lb,ub,n,m,cheby,suppl,xunits,yunits,lam,pi0);
% u = uv{1};
% Initialize random number seed to make the random number sequence equal to the case of global solution to the OCP
rng(42,'twister');
for i = 1:1200
    rand(1,2);
end
rand(2,1);
rand(2,1);
u = hessian_fmincon(nci,nphi,av,np,idx_x0,ti,ts0,tf0,x0,[],sc,lam,pi0,lb,ub,0.001); % Local solution to the OCP for the modified model
tsm0 = u(1); % Optimal switching time to arc with input at lower bound for the modified model
tfm0 = u(2); % Optimal final time for the modified model
end

% Draw plots (Figure 1 in the paper)
figure,hold on;
yyaxis left;
plot(1:niter,Jpk','-','LineWidth',1.5);
plot(1:niter,tfc*ones(1,niter)','--','LineWidth',1.5);
xlabel('$k$','Interpreter','LaTeX','FontSize',22);
ylabel('$\hat{\phi}^{p}(${\boldmath$\tau$}$_{k})$','Interpreter','LaTeX','FontSize',22);
set(gca,'FontSize',20,'XLim',[1,niter],'XTick',1:2:5);
yyaxis right;
plot(1:niter,Tpk','-','LineWidth',1.5);
plot(1:niter,zeros(1,niter)','--','LineWidth',1.5);
ylabel('$\hat{\psi_{j}}^{p}(${\boldmath$\tau$}$_{k})$','Interpreter','LaTeX','FontSize',22);
set(gca,'Position',[0.18,0.135,0.64,0.78]);
set(gca,'FontSize',20,'XLim',[1,niter],'XTick',1:2:5);
