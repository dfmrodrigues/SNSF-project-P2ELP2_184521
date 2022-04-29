% An implementation of MPC using direct methods
% Adapted from Joel Andersson, 2016

clc
close all
clear
rng(42)

import casadi.*

format long

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

% pm = pp;

% Initial solution using direct parsimonious input parameterization
pmp = pm;
Spmp = Spm;
av = [-2,1,-1];
nt = length(av);
ni = sum(av>0);
tf0 = tf;
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
deriv = 0;
lam0 = {};
[tf_optp,ts_optp,z0_optp,p0_optp,x0_optp,lam_optp,phi_optp,psi_optp,F_optp] =...
    direct_parsimonious(t0,tf,N,cvode,lbu,ubu,av,tf0,ts0,z00,p00,lbx,ubx,x0,t,x,u,p,pmp,xdot,qdot,phi,psi,eps0,deriv,lam0);
tv_optp = linspace(t0,tf_optp,N+1);
[J_optp,T_optp,x_optp] = calc_parsimonious(N,x0_optp(1:nx),tf_optp,ts_optp,z0_optp,p0_optp,pmp,phi_optp,psi_optp,F_optp);

% Simulate MPC
J_simp = zeros(1,N+1);
T_simp = zeros(npsi,N+1);
tf_simp = zeros(1,N+1);
ts_simp = zeros(nt-1,N+1);
z0_simp = zeros(ni*nu,N+1);
p0_simp = zeros(ni*nu,N+1);
tv_simp = zeros(1,N+1);
x_simp = zeros(size(x_optp));
tf0k = tf_optp;
ts0k = ts_optp;
z00k = z0_optp;
p00k = p0_optp;
x0k = x_optp(:,1);
t0k = t0;
lam0k = lam_optp;
for k = 1:(N+1)
    x_simp(:,k) = x0k(:,1);
    tv_simp(k) = t0k;
    while(true)
    [tf_optk,ts_optk,z0_optk,p0_optk,x0_optk,lam_optk,phi_optk,psi_optk,F_optk] =...
        direct_parsimonious(t0k,tf,N-k+1,cvode,lbu,ubu,av,tf0k,ts0k,z00k,p00k,lbx,ubx,x0k,t,x,u,p,pmp,xdot,qdot,phi,psi,eps0,deriv,lam0k);
    [J_optk,T_optk,x_optk] = calc_parsimonious(N-k+1,x0_optk(1:nx),tf_optk,ts_optk,x0_optk((nx+1):end),p0_optk,pmp,phi_optk,psi_optk,F_optk);
    if(any(T_optk>1e-3))
        ts0k(find(ts0k<t0k,1,'last')) = t0k+1e-3;
    else
        break;
    end
    end
    disp(sum(sum(abs(tf_optk-tf_optp))))
    disp(sum(sum(abs(ts_optk-ts_optp))))
    disp(sum(sum(abs(z0_optk-z0_optp))))
    disp(sum(sum(abs(p0_optk-p0_optp))))
    disp(sum(sum(abs(x_optk(:,1:end)-x_optp(:,k:end)))))
    J_simp(:,k) = J_optk;
    T_simp(:,k) = T_optk;
    tf_simp(:,k) = tf_optk;
    ts_simp(:,k) = ts_optk;
    z0_simp(:,k) = z0_optk;
    p0_simp(:,k) = p0_optk;
    if(k==N+1)
        break;
    end
    tf0k = tf_optk;
    ts0k = ts_optk;
    z00k = z0_optk;
    p00k = p0_optk;
    fk = F_optk{1}('x0',x_optk(:,1),'p',pp,'tf',tf_optk,'ts',ts_optk,'z0',z0_optk,'p0',p0_optk);
    x0k = full(fk.xf);
    t0k = t0k+(tf_optk-t0k)/(N-k+1);
    lam0k = lam_optk;
end
disp(tf_simp-tf_optp)
disp(ts_simp-ts_optp)
disp(z0_simp-z0_optp)
disp(p0_simp-p0_optp)
disp(x_simp-x_optp)

% Plot the solution
for c = 'os'
if(c=='o')
    tf_opt = tf_optp;
    ts_opt = ts_optp;
    z0_opt = z0_optp;
    p0_opt = p0_optp;
    tv_opt = tv_optp;
    x_opt = x_optp;
elseif(c=='s')
    tf_opt = tf_simp;
    ts_opt = ts_simp;
    z0_opt = z0_simp;
    p0_opt = p0_simp;
    tv_opt = tv_simp;
    x_opt = x_simp;
end
labels = cell(1,nx+nu);
figure;
hold on
for i = 1:nx
    plot(tv_opt,x_opt(i,:)-x0(i));
    labels{i} = ['x',num2str(i)];
end
for i = 1:nu
    labels{nx+i} = ['u',num2str(i)];
end
prev = [];
for k = 1:size(ts_opt,2)
if(c=='o')
ieb = 1:nt;
teb = [t0,ts_opt',tf_opt];
elseif(c=='s')
if(k==N+1)
    break;
end
teb = [t0,ts_opt(:,k)',tf_opt(k)];
ieb = (find(teb>tv_opt(k)+1e-3,1,'first')-1):find(teb<tv_opt(k+1)-1e-3,1,'last');
teb = [tv_opt(k),teb(ieb(2:end)),tv_opt(k+1)];
end
for i = 1:(length(teb)-1)
    if(av(ieb(i))>0)
        hold on;ax = gca;ax.ColorOrderIndex = nx+1;
        j = av(ieb(i));
        plot([teb(i),teb(i)],[prev,x_opt(nx+j,k)]);
        hold on;ax = gca;ax.ColorOrderIndex = nx+1;
        prev = x_opt(nx+j,k)+p0_opt(j,k)*(teb(i+1)-teb(i));
        plot([teb(i),teb(i+1)],[x_opt(nx+j,k),prev]);
    elseif(av(ieb(i))==-1)
        hold on;ax = gca;ax.ColorOrderIndex = nx+1;
        plot([teb(i),teb(i)],[prev,lbu(1)]);
        hold on;ax = gca;ax.ColorOrderIndex = nx+1;
        prev = lbu(1);
        plot([teb(i),teb(i+1)],[lbu(1),prev]);
    elseif(av(ieb(i))==-2)
        hold on;ax = gca;ax.ColorOrderIndex = nx+1;
        plot([teb(i),teb(i)],[prev,ubu(1)]);
        hold on;ax = gca;ax.ColorOrderIndex = nx+1;
        prev = ubu(1);
        plot([teb(i),teb(i+1)],[ubu(1),prev]);
    end
end
end
xlabel('t')
legend(labels)
end
keyboard

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
for c = 'os'
if(c=='o')
    tf_opt = tf_opts;
    u_opt = u_opts;
    tv_opt = tv_opts;
    x_opt = x_opts;
elseif(c=='s')
    tf_opt = tf_sims;
    u_opt = u_sims;
    tv_opt = tv_sims;
    x_opt = x_sims;
end
labels = cell(1,nx+nu);
figure;
hold on
for i = 1:nx
    plot(tv_opt,x_opt(i,:)-x0(i));
    labels{i} = ['x',num2str(i)];
end
for i = 1:nu
    stairs(tv_opt,u_opt(i,:))
    labels{nx+i} = ['u',num2str(i)];
end
xlabel('t')
legend(labels)
end

% Initial solution using direct multiple shooting
pmm = pm;
Spmm = Spm;
tf0 = tf;
u0 = (lbu+ubu)/2*ones(1,N);
x0 = x0(1:nx);
lam0 = {};
[tf_optm,u_optm,x0_optm,lam_optm,phi_optm,psi_optm,F_optm] =...
    direct_multiple_shooting(t0,tf,N,cvode,lbu,ubu,tf0,u0,lbx,ubx,x0,t,x,u,p,pmm,xdot,qdot,phi,psi,lam0);
tv_optm = linspace(t0,tf_optm,N+1);
[J_optm,T_optm,x_optm] = calc_multiple_shooting(N,x0_optm,tf_optm,u_optm,pmm,phi_optm,psi_optm,F_optm);

% Simulate MPC
J_simm = zeros(1,N+1);
T_simm = zeros(npsi,N+1);
tf_simm = zeros(size(tf_optm));
u_simm = zeros(size(u_optm));
tv_simm = zeros(1,N+1);
x_simm = zeros(size(x_optm));
tf0k = tf_optm;
u0k = u_optm;
x0k = x_optm(:,1);
t0k = t0;
lam0k = lam_optm;
for k = 1:(N+1)
    x_simm(:,k) = x0k(:,1);
    tv_simm(k) = t0k;
    [tf_optk,u_optk,x0_optk,lam_optk,phi_optk,psi_optk,F_optk] =...
        direct_multiple_shooting(t0k,tf,N-k+1,cvode,lbu,ubu,tf0k,u0k,lbx,ubx,x0k,t,x,u,p,pmm,xdot,qdot,phi,psi,lam0k);
    [J_optk,T_optk,x_optk] = calc_multiple_shooting(N-k+1,x0_optk,tf_optk,u_optk,pmm,phi_optk,psi_optk,F_optk);
    disp(sum(sum(abs(tf_optk-tf_optm))))
    disp(sum(sum(abs(u_optk(:,1:end-1)-u_optm(:,k:end-1)))))
    disp(sum(sum(abs(x_optk(:,1:end)-x_optm(:,k:end)))))
    J_simm(:,k) = J_optk;
    T_simm(:,k) = T_optk;
    tf_simm(:,k) = tf_optk;
    u_simm(:,k) = u_optk(:,1);
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
disp(tf_simm-tf_optm)
disp(u_simm-u_optm)
disp(x_simm-x_optm)

% Plot the solution
for c = 'os'
if(c=='o')
    tf_opt = tf_optm;
    u_opt = u_optm;
    tv_opt = tv_optm;
    x_opt = x_optm;
elseif(c=='s')
    tf_opt = tf_simm;
    u_opt = u_simm;
    tv_opt = tv_simm;
    x_opt = x_simm;
end
labels = cell(1,nx+nu);
figure;
hold on
for i = 1:nx
    plot(tv_opt,x_opt(i,:)-x0(i));
    labels{i} = ['x',num2str(i)];
end
for i = 1:nu
    stairs(tv_opt,u_opt(i,:))
    labels{nx+i} = ['u',num2str(i)];
end
xlabel('t')
legend(labels)
end
