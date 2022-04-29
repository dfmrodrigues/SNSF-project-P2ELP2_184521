% An implementation of direct single shooting
% Adapted from Joel Andersson, 2016

function [tf_opt,u_opt,x0_opt,lam_opt,phi_opt,psi_opt,F_opt] =...
    direct_single_shooting(t0,tf,N,cvode,lbu,ubu,tf0,u0,lbx,ubx,x0,t,x,u,p,pm,xdot,qdot,phi,psi,lam0)

import casadi.*

nx = size(x,1);
nu = size(u,1);
np = size(p,1);
npsi = size(psi,1);
if(~isempty(lam0))
    lam_tf0 = lam0{1};
    lam_u0 = lam0{2};
    lam_dt0 = lam0{3};
    lam_x0 = lam0{4};
    lam_dx0 = lam0{5};
    lam_T0 = lam0{6};
end
p = [p;u];
dt = SX.sym('dt');
p = [p;dt];

xdoti = dt*xdot;
qdoti = dt*qdot;
% Formulate discrete time dynamics
if cvode==true
   % CVODES from the SUNDIALS suite
   dae = struct('x',x,'p',p,'ode',xdoti,'quad',qdoti);
   opts = struct('tf',1,'reltol',1e-9,'abstol',1e-12);
   F = integrator('F','cvodes',dae,opts);
else
   % Fixed step Runge-Kutta 4 integrator
   M = 200; % RK4 steps per interval
   DT = 1/M;
   f = Function('f',{x,p},{xdoti,qdoti});
   X0 = MX.sym('X0',nx);
   P = MX.sym('P',np+nu+1);
   X = X0;
   Q = 0;
   for j = 1:M
       [k1_x,k1_q] = f(X,P);
       [k2_x,k2_q] = f(X+DT/2*k1_x,P);
       [k3_x,k3_q] = f(X+DT/2*k2_x,P);
       [k4_x,k4_q] = f(X+DT*k3_x,P);
       X = X+DT/6*(k1_x+2*k2_x+2*k3_x+k4_x);
       Q = Q+DT/6*(k1_q+2*k2_q+2*k3_q+k4_q);
    end
    F = Function('F',{X0,P},{X,Q},{'x0','p'},{'xf','qf'});
end

% Start with an empty NLP
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;
T = zeros(npsi,1);
g = {};
lbg = [];
ubg = [];
if(~isempty(lam0))
    lam_w0 = [];
    lam_g0 = [];
end

% Deal with free final time
Tf = MX.sym('Tf');
flagt = double(tf0>t0);
if(flagt)
    w = [w,{Tf}];
    lbw = [lbw;t0];
    ubw = [ubw;tf];
    w0 = [w0;tf0];
    if(~isempty(lam0))
        lam_w0 = [lam_w0;lam_tf0];
    end
else
    Tf = t0;
end
dT = (Tf-t0)/N;
g = [g,{-dT}];
lbg = [lbg;-inf];
ubg = [ubg;0];
if(~isempty(lam0))
    lam_g0 = [lam_g0;lam_dt0];
end

% Formulate the NLP
Xk = x0(:,1);
g = [g,{Xk}];
lbg = [lbg;-inf*ones(nx,1)];
ubg = [ubg;inf*ones(nx,1)];
for k = 1:N
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k-1)],nu);
    w = [w,{Uk}];
    lbw = [lbw;lbu];
    ubw = [ubw;ubu];
    w0 = [w0;u0(:,k)];

    % Integrate till the end of the interval
    Fk = F('x0',Xk,'p',[pm;Uk;dT]);
    Xk = Fk.xf;
    J = J+Fk.qf;

    % Add inequality constraint
    g = [g,{Xk}];
    lbg = [lbg;lbx];
    ubg = [ubg;ubx];
end
if(~isempty(lam0))
    lam_w0 = [lam_w0;reshape(lam_u0,[],1)];
    lam_w0 = lam_w0(1:(end-nu));
    lam_g0 = [lam_g0;lam_x0(:)+0*lam_dx0(:)];
end

% Add terminal constraints
phi_opt = Function('phi',{x,t},{phi},{'xf','tf'},{'phif'});
psi_opt = Function('psi',{x,t},{psi},{'xf','tf'},{'psif'});
Phif = phi_opt('xf',Xk,'tf',Tf);
J = J+Phif.phif;
Psif = psi_opt('xf',Xk,'tf',Tf);
T = T+Psif.psif;
g = [g,{T}];
lbg = [lbg;-inf*ones(npsi,1)];
ubg = [ubg;zeros(npsi,1)];
if(~isempty(lam0))
    lam_g0 = [lam_g0;lam_T0];
end

if(~isempty(w))
    % Create an NLP solver
    prob = struct('f',J,'x',vertcat(w{:}),'g',vertcat(g{:}));
    % opts = struct('error_on_fail',false,...
    %     'qpsol_options',struct('print_header',false,'print_iter',false,'print_info',false));
    % solver = nlpsol('solver','qrsqp',prob,opts);
    if(~isempty(lam0))
        opts = struct('ipopt',struct('print_level',5,'tol',1e-7,'dual_inf_tol',1e-7,...
            'warm_start_init_point','yes','mu_init',1e-9,'warm_start_bound_push',1e-12,...
            'warm_start_mult_bound_push',1e-12,'warm_start_slack_bound_push',1e-12));
    else
        opts = struct('ipopt',struct('print_level',5,'tol',1e-9,'dual_inf_tol',1e-9));
    end
    solver = nlpsol('solver','ipopt',prob,opts);
    
    % Solve the NLP
    if(~isempty(lam0))
        sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg,...
            'lam_x0',lam_w0,'lam_g0',lam_g0);
    else
        sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
    end
else
    sol = struct('f',J,'x',[],'g',vertcat(g{:}),'lam_x',[],'lam_g',lam_g0);
end

% Return the solution
w_opt = full(sol.x);
lam_w_opt = full(sol.lam_x);
lam_g_opt = full(sol.lam_g);
if(flagt)
    tf_opt = w_opt(1:flagt);
else
    tf_opt = tf0;
end
w_opt = reshape([w_opt((flagt+1):end);nan*ones(nu,1)],nu,[]);
u_opt = w_opt;
x0_opt = x0(:,1);
lam_tf_opt = lam_w_opt(1:flagt);
lam_w_opt = reshape([lam_w_opt((flagt+1):end);nan*ones(nu,1)],nu,[]);
lam_u_opt = lam_w_opt;
lam_dt_opt = lam_g_opt(1);
lam_x_opt = reshape(lam_g_opt(2:(end-npsi)),nx,[]);
lam_dx_opt = reshape(lam_g_opt(2:(end-npsi)),nx,[]);
lam_T_opt = lam_g_opt((end-npsi+1):end);
X0 = MX.sym('X0',nx);
P = MX.sym('P',np);
TF = MX.sym('TF',1);
U = MX.sym('U',nu);
F_opt = cell(1,N);
for k = 1:N
    tb = t0+(TF-t0)/N*(k-1);
    te = t0+(TF-t0)/N*k;
    Fk = F('x0',X0,'p',[P;U;te-tb]);
    Xk = Fk.xf;
    Qk = Fk.qf;
    F_opt{k} = Function(['F_',num2str(k)],{X0,P,TF,U},{Xk,Qk},{'x0','p','tf','u'},{'xf','qf'});
end
lam_opt = {lam_tf_opt,lam_u_opt,lam_dt_opt,lam_x_opt,lam_dx_opt,lam_T_opt};

end
