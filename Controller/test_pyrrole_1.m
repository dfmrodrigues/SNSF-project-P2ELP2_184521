clc
clear
close all
k10 = 0.053*10; % l mol-1 min-1
k20 = 0.128*10; % l mol-1 min-1
Ea1 = 20000; % kJ kmol-1
Ea2 = 10000; % kJ kmol-1
h = 0.02; % min
t_cl = 3; % min
q = 151;
deltat = (q-1)*h; %min
while(true)
    reply = input('Type "n" for fixed plant and design parameters, "k" for varying rate constants, "e" for varying activation energies, "t" for varying design parameters\n','s');
    if(strcmp(reply,'n')||strcmp(reply,'k')||strcmp(reply,'e')||strcmp(reply,'t'))
        break;
    end
end
if(strcmp(reply,'k'))
    k10 = 0.2:0.2:2;
    k20 = 0.4:0.4:4;
    X = ones(length(k20),1)*k10;
    Y = k20'*ones(1,length(k10));
    Z = zeros(length(k20),length(k10));
elseif(strcmp(reply,'e'))
    Ea1 = 10000:2000:30000;
    Ea2 = 5000:1000:15000;
    X = ones(length(Ea2),1)*Ea1;
    Y = Ea2'*ones(1,length(Ea1));
    Z = zeros(length(Ea2),length(Ea1));
elseif(strcmp(reply,'t'))
    t_cl = 1:0.5:10;
    q = 51:25:501;
    deltat = (q-1)*h;
    X = ones(length(deltat),1)*t_cl;
    Y = deltat'*ones(1,length(t_cl));
    Z = zeros(length(deltat),length(t_cl));
else
    Z = 0;
end
for j = 1:size(Z,1)
for k = 1:size(Z,2)
if(strcmp(reply,'k'))
    k10 = X(j,k);
    k20 = Y(j,k);
elseif(strcmp(reply,'e'))
    Ea1 = X(j,k);
    Ea2 = Y(j,k);
elseif(strcmp(reply,'t'))
    t_cl = X(j,k);
    deltat = Y(j,k);
end
Tref = 298.15; % K
R = 8.314; % kJ kmol-1 K-1
A(1) = k10*exp(Ea1/R/Tref); % l mol-1 min-1
A(2) = k20*exp(Ea2/R/Tref); % l mol-1 min-1
dH(1) = -70; % kJ mol-1
dH(2) = -100; % kJ mol-1
cin(1) = 2; % mol l-1
cin(2) = 1.5; % mol l-1
Tin(1) = 0; % kJ l-1
Tin(2) = 0; % kJ l-1
V = 500; % l
Rho = 0.87; % kg l-1
Cp = 155.96/92.14; % kJ kg-1 K-1
nCp = Rho*Cp*V; % kJ K-1
T0 = Tref;
Q0 = 0;
fA0 = 700/60; % l min-1
fB0 = 1100/60; % l min-1
s0 = solution_pyrrole_fA_fB(R,[Ea1,Ea2],A,dH,cin(1),cin(2),V,T0,fA0,fB0,0);
nA0 = s0(1);
nB0 = s0(2);
nC0 = s0(3);
nD0 = s0(4);
qex0 = s0(6);
r10 = k10*nA0*nB0/V;
r20 = k20*nB0^2/V;
n0 = [nA0;nB0;nC0;nD0];
u0 = [qex0;fA0;fB0];
d_ol = 0.005;
d_cl = 0.005;
d_u = d_ol*[qex0;fA0;fB0];
d_y = d_cl*[(T0-273.15)*nCp;nA0;nB0];
te = 90; % min

aQQ = -dH(1)*r10*(Ea1/R/T0^2/nCp)-dH(2)*r20*(Ea2/R/T0^2/nCp)-(fA0+fB0)/V;
aQnA = -dH(1)*k10*nB0/V;
aQnB = -dH(1)*k10*nA0/V-dH(2)*k20*2*nB0/V;
aQnC = 0;
aQnD = 0;
anAQ = -r10*(Ea1/R/T0^2/nCp);
anAnA = -(fA0+fB0)/V-k10*nB0/V;
anAnB = -k10*nA0/V;
anAnC = 0;
anAnD = 0;
anBQ = -r10*(Ea1/R/T0^2/nCp)-2*r20*(Ea2/R/T0^2/nCp);
anBnA = -k10*nB0/V;
anBnB = -(fA0+fB0)/V-k10*nA0/V-2*k20*2*nB0/V;
anBnC = 0;
anBnD = 0;
anCQ = r10*(Ea1/R/T0^2/nCp);
anCnA = k10*nB0/V;
anCnB = k10*nA0/V;
anCnC = -(fA0+fB0)/V;
anCnD = 0;
anDQ = r20*(Ea2/R/T0^2/nCp);
anDnA = 0;
anDnB = k20*2*nB0/V;
anDnC = 0;
anDnD = -(fA0+fB0)/V;
bQqex = 1;
bQfA = Tin(1)-Q0/V;
bQfB = Tin(2)-Q0/V;
bnAqex = 0;
bnAfA = cin(1)-nA0/V;
bnAfB = -nA0/V;
bnBqex = 0;
bnBfA = -nB0/V;
bnBfB = cin(2)-nB0/V;
bnCqex = 0;
bnCfA = -nC0/V;
bnCfB = -nC0/V;
bnDqex = 0;
bnDfA = -nD0/V;
bnDfB = -nD0/V;
A_ol = [aQQ,aQnA,aQnB,aQnC,aQnD;anAQ,anAnA,anAnB,anAnC,anAnD;anBQ,anBnA,anBnB,anBnC,anBnD;anCQ,anCnA,anCnB,anCnC,anCnD;anDQ,anDnA,anDnB,anDnC,anDnD];
B_ol = [bQqex,bQfA,bQfB;bnAqex,bnAfA,bnAfB;bnBqex,bnBfA,bnBfB;bnCqex,bnCfA,bnCfB;bnDqex,bnDfA,bnDfB];
C = eye(3,5);
S = eye(3);
u_cl = 1:3;
c_cl = 1:3;
n_cl = zeros(1,0);

n_u = size(B_ol,2);
n_x = size(A_ol,1);
n_y = size(C,1);
n_c = size(S,1);
dsady_ol = S*C*(-(fA0+fB0)/V*eye(size(A_ol)))*pinv(C);
dsadu_ol = S*C*B_ol;
Bc_ol = dsadu_ol;
C_p = C;
D_p = zeros(n_y,n_u);
C_ol = [S*C_p;zeros(n_u,n_x)];
D_ol = [S*D_p;eye(n_u)];
sys_ol = minreal(tf(ss(A_ol,B_ol,C_ol,D_ol)));
tau_ol = -1./real(eig(ssdata(ss(sys_ol,'min'))));

D_ay = dsady_ol;
D_au = dsadu_ol;
A_e = [-[6/deltat*eye(n_c);12/deltat^2*eye(n_c)],[eye(n_c);zeros(n_c)]];
B_ey = [12/deltat^2*S;zeros(n_c,n_y)];
B_ea = [zeros(n_c);-12/deltat^2*eye(n_c)];
C_e1 = [eye(n_c),zeros(n_c)];
n_z = 2*n_c;
D_ue1 = -Bc_ol\eye(n_u);
D_uy = -Bc_ol\(S/t_cl+dsady_ol);
D_ur0 = Bc_ol\eye(n_u)/t_cl;
A_k = A_e+B_ea*D_au*D_ue1*C_e1;
B_ky = B_ey+B_ea*D_ay+B_ea*D_au*D_uy;
B_kr0 = B_ea*D_au*D_ur0;
C_k = D_ue1*C_e1;
D_ky = D_uy;
D_kr0 = D_ur0;
A_cl = [A_ol+B_ol*D_ky*C_p,B_ol*C_k;B_ky*C_p,A_k];
B_clr0 = [B_ol*D_kr0;B_kr0];
B_cld0 = [B_ol;zeros(n_z,n_u)];
B_clw0 = [B_ol*D_ky;B_ky];
C_cl = [S*C_p,S*zeros(n_y,n_z);D_ky*C_p,C_k];
D_clr0 = [S*zeros(n_y,n_u);D_kr0];
D_cld0 = [S*zeros(n_y,n_u);eye(n_u)];
D_clw0 = [S*zeros(n_y,n_y);D_ky];
sys_cl = minreal(tf(ss(A_cl,eye(n_x+n_z,n_x),C_cl,[S*zeros(n_y,n_x);zeros(n_u,n_x)])));
tau_cl = -1./real(eig(ssdata(ss(sys_cl,'min'))));
sys_clr0 = minreal(tf(ss(A_cl,B_clr0,C_cl,D_clr0)));
tau_clr0 = -1./real(eig(ssdata(ss(sys_clr0,'min'))));
sys_cld0 = minreal(tf(ss(A_cl,B_cld0,C_cl,D_cld0)));
tau_cld0 = -1./real(eig(ssdata(ss(sys_cld0,'min'))));
sys_clw0 = minreal(tf(ss(A_cl,B_clw0,C_cl,D_clw0)));
tau_clw0 = -1./real(eig(ssdata(ss(sys_clw0,'min'))));
Z(j,k) = max(tau_clr0);
end
end
if(~strcmp(reply,'n'))
figure;
[C,ch] = contour(X,Y,Z);clabel(C,ch,'FontSize',18);
if(strcmp(reply,'k'))
set(gca,'FontSize',18,'XTick',0.2:0.2:2,'YTick',0.4:0.4:4);
title('$\tau_{cl}^{max}(\Delta t,\tau_c,${\boldmath$\theta$}$)$ [min]','Interpreter','LaTeX','FontSize',20)
xlabel('$k_{1,ref}$ [L mol$^{-1}$ min$^{-1}$]','Interpreter','LaTeX','FontSize',20);
ylabel('$k_{2,ref}$ [L mol$^{-1}$ min$^{-1}$]','Interpreter','LaTeX','FontSize',20);
elseif(strcmp(reply,'e'))
set(gca,'FontSize',18,'XTick',10000:2000:30000,'YTick',5000:1000:15000);
title('$\tau_{cl}^{max}(\Delta t,\tau_c,${\boldmath$\theta$}$)$ [min]','Interpreter','LaTeX','FontSize',20)
xlabel('$E_{a,1}$ [J mol$^{-1}$]','Interpreter','LaTeX','FontSize',20);
ylabel('$E_{a,2}$ [J mol$^{-1}$]','Interpreter','LaTeX','FontSize',20);
elseif(strcmp(reply,'t'))
set(gca,'FontSize',18,'XTick',1:1:10,'YTick',1:1:10);
title('$\tau_{cl}^{max}(\Delta t,\tau_c,${\boldmath$\theta$}$)$ [min]','Interpreter','LaTeX','FontSize',20)
xlabel('$\tau_c$ [min]','Interpreter','LaTeX','FontSize',20);
ylabel('$\Delta t$ [min]','Interpreter','LaTeX','FontSize',20);
end
else
Qr0 = @(i,t)Q0+(i==4)*d_cl*(T0-273.15)*nCp*ones(size(t));
nAr0 = @(i,t)nA0+(i==5)*d_cl*nA0*ones(size(t));
nBr0 = @(i,t)nB0+(i==6)*d_cl*nB0*ones(size(t));
vQ = @(i,t,Q)(Qr0(i,t)-Q)/t_cl;
vnA = @(i,t,n)(nAr0(i,t)-n(1))/t_cl;
vnB = @(i,t,n)(nBr0(i,t)-n(2))/t_cl;
qex = @(i,t,Q,n,yu1)(i<=3)*qex0+(i==1)*d_ol*qex0+...
    (i>3)*((-yu1(1)+vQ(i,t,Q))+(cin(2)*Q*(-yu1(2)+vnA(i,t,n))+cin(1)*Q*(-yu1(3)+vnB(i,t,n)))/(cin(1)*cin(2)*V-cin(2)*n(1)-cin(1)*n(2)));
fA = @(i,t,Q,n,yu1)(i<=3)*fA0+(i==2)*d_ol*fA0+...
    (i>3)*((cin(2)*V-n(2))*(-yu1(2)+vnA(i,t,n))+n(1)*(-yu1(3)+vnB(i,t,n)))/(cin(1)*cin(2)*V-cin(2)*n(1)-cin(1)*n(2));...
fB = @(i,t,Q,n,yu1)(i<=3)*fB0+(i==3)*d_ol*fB0+...
    (i>3)*(n(2)*(-yu1(2)+vnA(i,t,n))+(cin(1)*V-n(1))*(-yu1(3)+vnB(i,t,n)))/(cin(1)*cin(2)*V-cin(2)*n(1)-cin(1)*n(2));...
k = (1:(q-1))-1;
bk1 = [6*(q-k-1).*(k+1)/q/(q^2-1),0];
k = (1:q)-1;
ck1 = 12*(k-(q-1)/2)/q/(q^2-1)/h;
lags = deltat-k(1:(end-1))*h;
eu = @(i,t)[(i==7)*d_ol*qex0;(i==8)*d_ol*fA0;(i==9)*d_ol*fB0]*((t-[lags,0]>=0));
ey = @(i,t)[(i==10)*d_cl*(T0-273.15)*nCp;(i==11)*d_cl*nA0;(i==12)*d_cl*nB0;0;0]*((t-[lags,0]>=0).*(1-exp((-t+[lags,0])/6)));
ax = steps_1(sys_ol,sys_clr0,sys_cld0,sys_clw0,c_cl,n_cl,d_u,d_y,te);
for i = 1:12
    qexi = @(t,Q,n,yu1)qex(i,t,Q,n,yu1);
    fAi = @(t,Q,n,yu1)fA(i,t,Q,n,yu1);
    fBi = @(t,Q,n,yu1)fB(i,t,Q,n,yu1);
    eui = @(t)eu(i,t);
    eyi = @(t)ey(i,t);
    sol = dde23(@(t,y,Z)dynamics_pyrrole_1(t,y,Z,qexi,fAi,fBi,cin(1),cin(2),V,k10,k20,dH(1),dH(2),Ea1,Ea2,R,nCp,Tref,C,S,ck1,bk1,lags,eui,eyi,h,u0),...
        lags,[Q0;nA0;nB0;nC0;nD0;zeros(3,1)],0:h:te,ddeset('AbsTol',1e-9,'RelTol',1e-6));
    tout = 0:h:te;
    yout = interp1(sol.x,sol.y',tout)';
    ypout = interp1(sol.x,sol.yp',tout)';
    Q = yout(1,:);
    n = yout(2:5,:);
    u = ypout(6:8,:);
    figure(i);
    hold(ax{1,i},'on'),plot(ax{1,i},tout,Q-Q0,'g','LineWidth',2);set(ax{1,i},'Children',flipud(get(ax{1,i},'Children')));
    hold(ax{1,i},'on'),plot(ax{1,i},tout,Qr0(i,tout)-Q0,'r');
    title(ax{1,i},'');ylabel(ax{1,i},'$\delta Q(t)$ [kJ]','Interpreter','LaTeX','FontSize',20);
    p = xlabel(ax{1,i},'');set(p,'String','','Interpreter','LaTeX','FontSize',20);
    set(ax{1,i},'FontSize',18,'XTick',0:30:90,'XTickLabel',[]);
    hold(ax{2,i},'on'),plot(ax{2,i},tout,n(1,:)-nA0,'g','LineWidth',2);set(ax{2,i},'Children',flipud(get(ax{2,i},'Children')));
    hold(ax{2,i},'on'),plot(ax{2,i},tout,nAr0(i,tout)-nA0,'r');
    title(ax{2,i},'');ylabel(ax{2,i},'$\delta n_{\rm A}(t)$ [mol]','Interpreter','LaTeX','FontSize',20);
    p = xlabel(ax{2,i},'');set(p,'String','','Interpreter','LaTeX','FontSize',20);
    set(ax{2,i},'FontSize',18,'XTick',0:30:90,'XTickLabel',[]);
    hold(ax{3,i},'on'),plot(ax{3,i},tout,n(2,:)-nB0,'g','LineWidth',2);set(ax{3,i},'Children',flipud(get(ax{3,i},'Children')));
    hold(ax{3,i},'on'),plot(ax{3,i},tout,nBr0(i,tout)-nB0,'r');
    title(ax{3,i},'');ylabel(ax{3,i},'$\delta n_{ \rm B}(t)$ [mol]','Interpreter','LaTeX','FontSize',20);
    p = xlabel(ax{3,i},'');set(p,'String','','Interpreter','LaTeX','FontSize',20);
    set(ax{3,i},'FontSize',18,'XTick',0:30:90,'XTickLabel',[]);
    hold(ax{4,i},'on'),plot(ax{4,i},tout,u(1,:)-qex0,'g','LineWidth',2);set(ax{4,i},'Children',flipud(get(ax{4,i},'Children')));
    title(ax{4,i},'');ylabel(ax{4,i},'$\delta q_{ex}(t)$ [kJ min$^{-1}$]','Interpreter','LaTeX','FontSize',20);
    p = xlabel(ax{4,i},'');set(p,'String','','Interpreter','LaTeX','FontSize',20);
    set(ax{4,i},'FontSize',18,'XTick',0:30:90,'XTickLabel',[]);
    hold(ax{5,i},'on'),plot(ax{5,i},tout,u(2,:)-fA0,'g','LineWidth',2);set(ax{5,i},'Children',flipud(get(ax{5,i},'Children')));
    title(ax{5,i},'');ylabel(ax{5,i},'$\delta F_{\rm A}(t)$ [L min$^{-1}$]','Interpreter','LaTeX','FontSize',20);
    p = xlabel(ax{5,i},'');set(p,'String','$t$ [min]','Interpreter','LaTeX','FontSize',20);
    set(ax{5,i},'FontSize',18,'XTick',0:30:90,'XTickLabel',0:30:90);
    hold(ax{6,i},'on'),plot(ax{6,i},tout,u(3,:)-fB0,'g','LineWidth',2);set(ax{6,i},'Children',flipud(get(ax{6,i},'Children')));
    title(ax{6,i},'');ylabel(ax{6,i},'$\delta F_{\rm B}(t)$ [L min$^{-1}$]','Interpreter','LaTeX','FontSize',20);
    p = xlabel(ax{6,i},'');set(p,'String','$t$ [min]','Interpreter','LaTeX','FontSize',20);
    set(ax{6,i},'FontSize',18,'XTick',0:30:90,'XTickLabel',0:30:90);
    set(gcf,'Units','normalized','OuterPosition',[0.1,0.1,0.8,0.8]);
end
end