function [Tf,CEMf] = appj_track(ts,tf)
% Define constants for the control approach
h = 0.02; % s
t_cl = 1; % s
q = 151;
deltat = (q-1)*h; % s
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
T0 = 310.15; % K
Pmin = 1.0; % W
Pmax = 5.0; % W
Tmax = 316.15; % K
CEMsp = 1.5; % min
k1 = m2*eta/Rho/Cp/pi/r^2/d*0.8;
te = 120; % s
C = 1;
S = 1;
% Specify control laws
Tr0 = @(t)Tmax*ones(size(t)); % Temperature reference for the constraint-seeking arc
vT = @(t,T)(Tr0(t)-T)/t_cl; % Desired temperature variation for the constraint-seeking arc
Ba = @(T)k1*ones(size(T)); % Value of Ba for the control approach
P = @(t,T,flagT,yu1)(flagT<=0.1)*Pmax+(flagT>0.1&&t<=ts)*(-yu1(1)+vT(t,T))/k1+(t>ts)*Pmin; % Control law for the manipulated variable \tilde{P} (Pmax, then constraint-seeking arc, then Pmin)
k = (1:(q-1))-1;
bk1 = [6*(q-k-1).*(k+1)/q/(q^2-1),0]; % FIR filter weights for y_c in the control approach
k = (1:q)-1;
ck1 = 12*(k-(q-1)/2)/q/(q^2-1)/h; % FIR filter weights for y_a in the control approach
lags = deltat-k(1:(end-1))*h; % Needed to implement FIR filter via the dde23f function below
% Specify disturbances
noiseu = normrnd(zeros(size(0:0.2:te)),0.1);
noisey = normrnd(zeros(size(0:0.2:te)),0.1);
eu = @(t)binterp_mex(0:0.2:te,noiseu,t); % Input disturbance
ey = @(t)binterp_mex(0:0.2:te,noisey,t); % Output disturbance
% Simulate the true APPJ subject to disturbances with the specified control laws
sol = dde23f(@(t,y,Z)dynamics_appj(t,y,Z,P,Ba,Tinf,Tb,Tmax,Rho,Cp,d,r,eta,K,m1,m2,C,S,ck1,bk1,lags,eu,ey,h,T0),...
    lags,[0;T0;0;0;0],0:h:te,ddeset('AbsTol',1e-9,'RelTol',1e-6)); % dde23f is a faster implementation of dde23
% Obtain final temperature and CEM and data for plots
tout = [0:h:tf,tf];
yout = interp1(sol.x,sol.y',tout)';
ypout = interp1(sol.x,sol.yp',tout)';
CEM = yout(1,:);
T = yout(2,:);
tP = ypout(4,:);
Tf = T(end)-T0;
CEMf = CEMsp-CEM(end);
% Draw plots (Figure 2 in the paper)
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
hold(ax{3},'on'),plot(ax{3},tout,tP(1,:),'b','LineWidth',2);set(ax{3},'Children',flipud(get(ax{3},'Children')));
hold(ax{3},'on'),plot(ax{3},tout,Pmin*ones(size(tout)),'r');
hold(ax{3},'on'),plot(ax{3},tout,Pmax*ones(size(tout)),'r');
title(ax{3},'');ylabel(ax{3},'$\tilde{P}(t)$ [W]','Interpreter','LaTeX','FontSize',20);
p = xlabel(ax{3},'');set(p,'String','$t$ [s]','Interpreter','LaTeX','FontSize',20);
set(ax{3},'FontSize',18,'XTick',0:30:90,'XTickLabel',0:30:90);
set(gcf,'Units','normalized'),set(gcf,'OuterPosition',[1,1,2,1].*get(gcf,'OuterPosition'));
end