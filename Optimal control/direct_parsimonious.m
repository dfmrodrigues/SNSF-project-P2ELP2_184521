% An implementation of direct parsimonious input parameterization
% Adapted from Joel Andersson, 2016

function [tf_opt,ts_opt,z0_opt,p0_opt,x0_opt,lam_opt,phi_opt,psi_opt,F_opt] =...
    direct_parsimonious(t0,tf,N,cvode,lbu,ubu,av,tf0,ts0,z00,p00,lbx,ubx,x0,t,x,u,p,pm,xdot,qdot,phi,psi,eps0,deriv,lam0)

import casadi.*

nx = size(x,1);
nu = size(u,1);
np = size(p,1);
npsi = size(psi,1);
if(length(eps0)==5)
    eps_J0 = eps0{1};
    eps_T0 = eps0{2};
    eps_dJdu0 = eps0{3};
    eps_dTdu0 = eps0{4};
    u0 = eps0{5};
else
    y = eps0{1};
    S = eps0{2};
    cm = eps0{3};
    lambdam = eps0{4};
    u0 = eps0{5};
    Uk = eps0{6};
    sc = eps0{7};
end
if(~isempty(lam0))
    lam_tf0 = lam0{1};
    lam_ts0 = lam0{2};
    lam_z00 = lam0{3};
    lam_p00 = lam0{4};
    lam_dt0 = lam0{5};
    lam_zf0 = lam0{6};
    lam_xi0 = lam0{7};
    lam_T0 = lam0{8};
end
nt = length(av);
ni = sum(av>0);
it0 = 1;
ii0 = 1+(t0>0);
while(it0<nt&&(tf0<=t0||t0>ts0(it0)))
    it0 = it0+1;
    ii0 = ii0+1;
end
nt0 = length(av(it0:end));
if(deriv>0)
    nzi0 = sum(av(it0:end)>0);
elseif(deriv==0)
    nzi0 = sum(av(ii0:end)>0);
end
npi0 = sum(av(ii0:end)>0);
ci = [1:ni,-(1:(length(lbu)+length(ubu)))];
crhs = cell(1,ni+length(lbu)+length(ubu));
prhs = cell(1,ni+length(lbu)+length(ubu));
dt = SX.sym('dt');
p = [p;dt];
for i = 1:ni
    s = [];
    for k = 1:nu
        xi = SX.sym(['x',num2str(nx+(i-1)*nu+k)]);
        s = [s;xi];
        x = [x;xi];
        lbx(nx+(i-1)*nu+k) = -inf;
        ubx(nx+(i-1)*nu+k) = inf;
    end
    crhs{i} = s;
    s = [];
    for j = 1:ni
        for k = 1:nu
            if(i==j)
                p0i = SX.sym(['p',num2str((i-1)*nu+k)]);
                s = [s;p0i];
                p = [p;p0i];
            else
                s = [s;0];
            end
        end
    end
    prhs{i} = s;
end
for i = 1:length(lbu)
    s = zeros(nu,1);
    for j = 1:nu
        if(i==j)
            s(j) = lbu(i);
        else
            s(j) = 0;
        end
    end
    crhs{ni+i} = s;
    prhs{ni+i} = zeros(ni*nu,1);
end
for i = 1:length(ubu)
    s = zeros(nu,1);
    for j = 1:nu
        if(i==j)
            s(j) = ubu(i);
        else
            s(j) = 0;
        end
    end
    crhs{ni+length(lbu)+i} = s;
    prhs{ni+length(lbu)+i} = zeros(ni*nu,1);
end

F0 = cell(1,nt);
for i = it0:nt
cases = find(av(i)==ci);
xdoti = substitute(dt*[xdot;prhs{cases}],u(:),crhs{cases});
qdoti = substitute(dt*qdot,u(:),crhs{cases});
% Formulate discrete time dynamics
if cvode==true
   % CVODES from the SUNDIALS suite
   dae = struct('x',x,'p',p,'ode',xdoti,'quad',qdoti);
   opts = struct('tf',1,'reltol',1e-9,'abstol',1e-12);
   F = integrator('F','cvodes',dae,opts);
else
   % Fixed step Runge-Kutta 4 integrator
   M = 500; % RK4 steps per interval
   DT = 1/M;
   f = Function('f',{x,p},{xdoti,qdoti});
   X0 = MX.sym('X0',nx+ni*nu);
   P = MX.sym('P',np+1+ni*nu);
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
F0{i} = F;
end

% Start with an empty NLP
w = {};
w0 = [];
lbw = [];
ubw = [];
if(length(eps0)==5)
    J = eps_J0;
    T = eps_T0;
else
    sumsqs = 0;
    sumsqv = zeros(size(Uk,2),1);
end
g = {};
lbg = [];
ubg = [];
if(~isempty(lam0))
    lam_w0 = [];
    lam_g0 = [];
end

% Deal with free final and switching times
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
dT = cell(1,nt);
for i = 1:nt
    if(i>it0)
        Tprev = Ti;
    else
        Tprev = t0;
    end
    if(i<nt&&i>=it0)
        Ti = MX.sym(['T_' num2str(i)]);
        if(length(eps0)==5)
            J = J+eps_dJdu0(:,i)*(Ti-u0(i));
            T = T+eps_dTdu0(:,i)*(Ti-u0(i));
        else
            sumsqv = sumsqv+(((Ti-u0(i)-Uk(i+1,:))/sc(i)).^2)';
        end
        w = [w,{Ti}];
        lbw = [lbw;-inf];
        ubw = [ubw;inf];
        w0 = [w0;ts0(i)];
        if(~isempty(lam0))
            lam_w0 = [lam_w0;lam_ts0(i)];
        end
    elseif(i<it0)
        Ti = t0;
    else
        Ti = Tf;
    end
    dT{i} = Ti-Tprev;
    g = [g,{-dT{i}}];
    lbg = [lbg;-inf];
    if(i<it0)
        ubg = [ubg;inf];
    elseif(av(i)>0&&t0==0&&deriv>0)
        ubg = [ubg;-min(ubu-lbu)/deriv];
    else
        ubg = [ubg;0];
    end
    if(~isempty(lam0))
        lam_g0 = [lam_g0;lam_dt0(end-nt+i)];
    end
end

% Formulate the NLP
Xi = x0(1:nx);
U00 = cell(1,nt);
for i = 1:nt
    if(av(i)>0)
        idx_i = ((av(i)-1)*nu+1):(av(i)*nu);
        if((deriv>0&&i>=it0)||(deriv==0&&i>=ii0))
            U0i = MX.sym(['U0_' num2str(i)],nu);
            if(length(eps0)==5)
                J = J+eps_dJdu0(:,nt-1+idx_i)*(U0i-u0(nt-1+idx_i));
                T = T+eps_dTdu0(:,nt-1+idx_i)*(U0i-u0(nt-1+idx_i));
            else
                sumsqv = sumsqv+(((U0i-u0(nt-1+idx_i)-Uk(nt-1+idx_i+1,:))/sc(nt-1+idx_i)).^2)';
            end
            w = [w,{U0i}];
            lbw = [lbw;lbu];
            ubw = [ubw;ubu];
            w0 = [w0;z00(idx_i)];
            if(~isempty(lam0))
                lam_w0 = [lam_w0;lam_z00(idx_i)];
            end
            U00{i} = U0i;
        else
            U0i = x0(nx+idx_i);
        end
        Xi = [Xi;U0i];
    end
end

P0 = [];
for i = 1:nt
    if(av(i)>0)
        idx_i = ((av(i)-1)*nu+1):(av(i)*nu);
        if(i>=ii0)
            Pi = MX.sym(['P_' num2str(i)],nu);
            if(length(eps0)==5)
                J = J+eps_dJdu0(:,nt-1+ni*nu+idx_i)*(Pi-u0(nt-1+ni*nu+idx_i));
                T = T+eps_dTdu0(:,nt-1+ni*nu+idx_i)*(Pi-u0(nt-1+ni*nu+idx_i));
            else
                sumsqv = sumsqv+(((Pi-u0(nt-1+ni*nu+idx_i)-Uk(nt-1+ni*nu+idx_i+1,:))/sc(nt-1+ni*nu+idx_i)).^2)';
            end
            w = [w,{Pi}];
            lbw = [lbw;-inf*ones(nu,1)];
            ubw = [ubw;inf*ones(nu,1)];
            w0 = [w0;p00(idx_i)];
            if(~isempty(lam0))
                lam_w0 = [lam_w0;lam_p00(idx_i)];
            end
            if(deriv>0)
                g = [g,{U00{i}+Pi*dT{i}}];
                lbg = [lbg;lbu];
                ubg = [ubg;ubu];
            elseif(deriv==0)
                g = [g,{Pi}];
                lbg = [lbg;zeros(nu,1)];
                ubg = [ubg;zeros(nu,1)];
            end
            if(~isempty(lam0))
                lam_g0 = [lam_g0;lam_zf0(end-ni*nu+idx_i)];
            end
        else
            Pi = p00(idx_i);
        end
        P0 = [P0;Pi];
    end
end
if(length(eps0)~=5)
    for l = 1:(npsi+1)
        Ks{l} = cm{l}*exp(-sumsqs/2/lambdam{l}^2);
        Kv{l} = cm{l}*exp(-sumsqv/2/lambdam{l}^2);
    end
    J = Kv{1}'/S{1}*y{1}+sqrt(Ks{1}-Kv{1}'/S{1}*Kv{1});
    T = [];
    for l = 1:npsi
        T = [T;Kv{l+1}'/S{l+1}*y{l+1}+sqrt(Ks{l+1}-Kv{l+1}'/S{l+1}*Kv{l+1})];
    end
end

for i = 1:nt
    if(i>=it0)
        % Integrate till the end of the interval
        Fi = F0{i}('x0',Xi,'p',[pm;dT{i};P0]);
        Xi = Fi.xf;
        J = J+Fi.qf;
        
        g = [g,{Xi}];
        if(i>it0||av(i)<=0||(all(x0>=lbx)&&all(x0<=ubx)))
            lbg = [lbg;lbx];
            ubg = [ubg;ubx];
        else
            lbg = [lbg;-inf*ones(nx+ni*nu,1)];
            ubg = [ubg;inf*ones(nx+ni*nu,1)];
        end
    else
        g = [g,{zeros(nx+ni*nu,1)}];
        lbg = [lbg;-inf*ones(nx+ni*nu,1)];
        ubg = [ubg;inf*ones(nx+ni*nu,1)];
    end
    if(~isempty(lam0))
        lam_g0 = [lam_g0;lam_xi0(end+(i-1-nt)*(nx+ni*nu)+(1:(nx+ni*nu)))];
    end
end

% Add terminal constraints
phi_opt = Function('phi',{x,t},{phi},{'xf','tf'},{'phif'});
psi_opt = Function('psi',{x,t},{psi},{'xf','tf'},{'psif'});
Phif = phi_opt('xf',Xi,'tf',Tf);
J = J+Phif.phif;
Psif = psi_opt('xf',Xi,'tf',Tf);
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
        opts = struct('ipopt',struct('print_level',5,'tol',1e-3,'dual_inf_tol',1e-3,...
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
lam_tf_opt = lam_w_opt(1:flagt);
if(~isempty(lam0))
    ts_opt = ts0(:);
    z0_opt = z00(:);
    p0_opt = p00(:);
    lam_ts_opt = lam_ts0;
    lam_z0_opt = lam_z00;
    lam_p0_opt = lam_p00;
else
    ts_opt = zeros(nt-1,1);
    z0_opt = zeros(ni*nu,1);
    p0_opt = zeros(ni*nu,1);
    lam_ts_opt = zeros(nt-1,1);
    lam_z0_opt = zeros(ni*nu,1);
    lam_p0_opt = zeros(ni*nu,1);
end
x0_opt = x0;
for i = it0:(nt-1)
    ts_opt(i) = w_opt(flagt+i-nt+nt0);
    lam_ts_opt(i) = lam_w_opt(flagt+i-nt+nt0);
end
for i = 1:(nt-1)
    if(av(i)>0)
        idx_i = ((av(i)-1)*nu+1):(av(i)*nu);
        if((deriv>0&&i>=it0)||(deriv==0&&i>=ii0))
            z0_opt(idx_i) = w_opt(flagt+nt0-1+av(i)-ni+nzi0);
            lam_z0_opt(idx_i) = lam_w_opt(flagt+nt0-1+av(i)-ni+nzi0);
            x0_opt(nx+idx_i) = z0_opt(idx_i);
        end
        if(i>=ii0)
            p0_opt(idx_i) = w_opt(flagt+nt0-1+nzi0*nu+av(i)-ni+npi0);
            lam_p0_opt(idx_i) = lam_w_opt(flagt+nt0-1+nzi0*nu+av(i)-ni+npi0);
        end
    end
end
lam_dt_opt = lam_g_opt(1:nt);
lam_zf_opt = lam_g_opt((nt+1):(end-npsi-(nx+ni*nu)*nt));
lam_xi_opt = lam_g_opt((end-npsi-(nx+ni*nu)*nt+1):(end-npsi));
lam_T_opt = lam_g_opt((end-npsi+1):end);
X0 = MX.sym('X0',nx+ni*nu);
P = MX.sym('P',np);
TF = MX.sym('TF',1);
TS = MX.sym('TS',nt-1);
Z0 = MX.sym('Z0',ni*nu);
P0 = MX.sym('P0',ni*nu);
F_opt = cell(1,N);
i = it0;
for k = 1:N
    Xk = X0;
    Qk = 0;
    tb = t0+(TF-t0)/N*(k-1);
    te = t0+(TF-t0)/N*k;
    while(i<nt&&t0+(tf_opt-t0)/N*k>ts_opt(i))
        Fk = F0{i}('x0',Xk,'p',[P;TS(i)-tb;P0]);
        Xk = Fk.xf;
        Qk = Qk+Fk.qf;
        tb = TS(i);
        i = i+1;
    end
    Fk = F0{i}('x0',Xk,'p',[P;te-tb;P0]);
    Xk = Fk.xf;
    Qk = Qk+Fk.qf;
    F_opt{k} = Function(['F_',num2str(k)],{X0,P,TF,TS,Z0,P0},{Xk,Qk},{'x0','p','tf','ts','z0','p0'},{'xf','qf'});
end
lam_opt = {lam_tf_opt,lam_ts_opt,lam_z0_opt,lam_p0_opt,lam_dt_opt,lam_zf_opt,lam_xi_opt,lam_T_opt};

end
