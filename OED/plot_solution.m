%% Plot final Solution
%Uses the solution
clearvars
format long 
close all
clc

%% Define some nice colors for plotting
custom_cols = [
    0.7294    0.1412    0.0510;...
    0.0824    0.5216    0.0157;...
    0.1882    0.0902    0.8196;...
    0.9020    0.5294    0.1098;...
    0.6078    0.1882    0.8902];
%custom_cols(1,:) ~ deepred
%custom_cols(2,:) ~ deepgreen
%custom_cols(3,:) ~ deepblue
%custom_cols(4,:) ~ deeporange
%custom_cols(5,:) ~ deeppurple

%% Define design parameters
arcs = 3; %Arc sequence contains 3 arcs 
ncoeff = 2; %General case is where each arc has up to 2 coeffs (linear)

Dtrain = zeros(arcs,ncoeff+1); %Matrix with design vars

while(true)
    reply_Dtrain = input('Type the arc sequence:\n');
    if(all(reply_Dtrain==0)||isempty(reply_Dtrain)||all(reply_Dtrain==[1,-1,-2])||all(reply_Dtrain==[1,-2,-1])||...
            all(reply_Dtrain==[-1,1,-1])||all(reply_Dtrain==[-1,1,-2])||all(reply_Dtrain==[-2,1,-1])||all(reply_Dtrain==[-2,1,-2]))
        break;
    end
end
while(true)
    reply_betad = input(['Type b for beta distribution, ',...
        'type n for normal distribution:\n'],'s');
    if(strcmp(reply_betad,'b')||strcmp(reply_betad,'n'))
        break;
    end
end
betad = strcmp(reply_betad,'b');
while(true)
    reply_fast = input(['Type f for fast and coarse approximations of the expected utility, ',...
        'type s for slow and precise approximations of the expected utility:\n'],'s');
    if(strcmp(reply_fast,'f')||strcmp(reply_fast,'s'))
        break;
    end
end
fast = strcmp(reply_fast,'f');
while(true)
    reply_chance = input(['Type y for chance constraints, ',...
        'type n for no chance constraints:\n'],'s');
    if(strcmp(reply_chance,'y')||strcmp(reply_chance,'n'))
        break;
    end
end
chance = strcmp(reply_chance,'y');


%% Optimal Solution
if(all(reply_Dtrain==0))
    Dtrain = zeros(31,ncoeff+1);
    Dtrain(:,1) = (idinput(31,'prbs')+1)/2;
    Dtrain(:,ncoeff+1) = linspace(0,12-12/31,31);
elseif(~betad)
if(isempty(reply_Dtrain))
    Dtrain(1,1) = 0.000000000000000;
    Dtrain(2,1) = 0.736253127427679;
    Dtrain(2,2) = -0.165755283328793;
    Dtrain(3,1) = 1.000000000000000;
    Dtrain(2,ncoeff+1) = 4.186312557738331;
    Dtrain(3,ncoeff+1) = 7.367536075801477;
elseif(all(reply_Dtrain==[1,-1,-2]))
    Dtrain(1,1) = 1.0000;
    Dtrain(1,2) = 0.0000;
    Dtrain(2,1) = 0.0000;
    Dtrain(3,1) = 1.0000;
    Dtrain(2,ncoeff+1) = 2.3197;
    Dtrain(3,ncoeff+1) = 12.0000;
elseif(all(reply_Dtrain==[1,-2,-1]))
    Dtrain(1,1) = 0.4817;
    Dtrain(1,2) = -0.0903;
    Dtrain(2,1) = 1.0000;
    Dtrain(3,1) = 0.0000;
    Dtrain(2,ncoeff+1) = 5.3340;
    Dtrain(3,ncoeff+1) = 9.4772;
elseif(all(reply_Dtrain==[-1,1,-1]))
    Dtrain(1,1) = 0.0000;
    Dtrain(2,1) = 1.0000;
    Dtrain(2,2) = 0.0000;
    Dtrain(3,1) = 0.0000;
    Dtrain(2,ncoeff+1) = 5.1300;
    Dtrain(3,ncoeff+1) = 10.1582;
elseif(all(reply_Dtrain==[-1,1,-2]))
    Dtrain(1,1) = 0.0000;
    Dtrain(2,1) = 1.0000;
    Dtrain(2,2) = -0.1291;
    Dtrain(3,1) = 1.0000;
    Dtrain(2,ncoeff+1) = 4.9775;
    Dtrain(3,ncoeff+1) = 12.0000;
elseif(all(reply_Dtrain==[-2,1,-1]))
    Dtrain(1,1) = 1.0000;
    Dtrain(2,1) = 0.0000;
    Dtrain(2,2) = 0.1339;
    Dtrain(3,1) = 0.0000;
    Dtrain(2,ncoeff+1) = 1.7936;
    Dtrain(3,ncoeff+1) = 9.2595;
elseif(all(reply_Dtrain==[-2,1,-2]))
    Dtrain(1,1) = 1.0000;
    Dtrain(2,1) = 0.0000;
    Dtrain(2,2) = 0.0000;
    Dtrain(3,1) = 1.0000;
    Dtrain(2,ncoeff+1) = 2.3197;
    Dtrain(3,ncoeff+1) = 12.0000;
end
else
if(isempty(reply_Dtrain))
    Dtrain(1,1) = 0.000000000000000;
    Dtrain(2,1) = 0.736253127427679;
    Dtrain(2,2) = -0.165755283328793;
    Dtrain(3,1) = 1.000000000000000;
    Dtrain(2,ncoeff+1) = 4.186312557738331;
    Dtrain(3,ncoeff+1) = 7.367536075801477;
elseif(all(reply_Dtrain==[1,-1,-2]))
    Dtrain(1,1) = 1.0000;
    Dtrain(1,2) = 0.0000;
    Dtrain(2,1) = 0.0000;
    Dtrain(3,1) = 1.0000;
    Dtrain(2,ncoeff+1) = 2.2512;
    Dtrain(3,ncoeff+1) = 12.0000;
elseif(all(reply_Dtrain==[1,-2,-1]))
    Dtrain(1,1) = 0.5378;
    Dtrain(1,2) = -0.1029;
    Dtrain(2,1) = 1.0000;
    Dtrain(3,1) = 0.0000;
    Dtrain(2,ncoeff+1) = 5.2057;
    Dtrain(3,ncoeff+1) = 9.3303;
elseif(all(reply_Dtrain==[-1,1,-1]))
    Dtrain(1,1) = 0.0000;
    Dtrain(2,1) = 1.0000;
    Dtrain(2,2) = 0.0000;
    Dtrain(3,1) = 0.0000;
    Dtrain(2,ncoeff+1) = 5.1155;
    Dtrain(3,ncoeff+1) = 10.3411;
elseif(all(reply_Dtrain==[-1,1,-2]))
    Dtrain(1,1) = 0.0000;
    Dtrain(2,1) = 1.0000;
    Dtrain(2,2) = -0.0514;
    Dtrain(3,1) = 1.0000;
    Dtrain(2,ncoeff+1) = 5.3782;
    Dtrain(3,ncoeff+1) = 12.0000;
elseif(all(reply_Dtrain==[-2,1,-1]))
    Dtrain(1,1) = 1.0000;
    Dtrain(2,1) = 0.1044;
    Dtrain(2,2) = 0.1153;
    Dtrain(3,1) = 0.0000;
    Dtrain(2,ncoeff+1) = 1.5561;
    Dtrain(3,ncoeff+1) = 9.3240;
elseif(all(reply_Dtrain==[-2,1,-2]))
    Dtrain(1,1) = 1.0000;
    Dtrain(2,1) = 0.2359;
    Dtrain(2,2) = -0.0240;
    Dtrain(3,1) = 1.0000;
    Dtrain(2,ncoeff+1) = 2.1886;
    Dtrain(3,ncoeff+1) = 12.0000;
end
end

k = 0.2;
K = k^2*eye(2);
if(~betad)
if(~fast)
w = [
    0.0201051368113866*ones(30,1)
    0.0119139564145941*ones(30,1)
    0.0012962636171673*ones(30,1)
    0.0000179764901854*ones(30,1)];% Normally distributed prior, 120 points
theta = [
    0.160635659707143*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]';...
    0.373712306584429*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]';...
    0.602436406397929*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]';...
    0.866951943872305*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]'];% Normally distributed prior, 120 points
else
w = [
   0.210491910111102 0.210491910111102 0.210491910111102 0.210491910111102 0.027777777777778 0.027777777777778...
   0.027777777777778 0.027777777777778 0.011730312111120 0.011730312111120 0.011730312111120 0.011730312111120]';% Normally distributed prior, 12 points
theta = [
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
   0.396335765891742   0.396335765891742];% Normally distributed prior, 12 points
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
theta = [th1(:),th2(:)];
wf = @(k,thetak)(2*np+a(k)+b(k)-1)./((c(k)^2-thetak.^2).*double(subs(diff(orth_pols(k,np),sym('x')),sym('x'),thetak)).^2);
w1gv = wf(1,th1gv);
w2gv = wf(2,th2gv);
[w1,w2] = meshgrid(w1gv,w2gv);
w = w1(:).*w2(:);
else
w = [
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
   0.016188049850202];% Beta distributed prior, 19 points
theta = [
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
   0.440109933805545   0.210778895104562];% Beta distributed prior, 19 points
end
end
sigma_z = 0.1;
m_z = 10000;
Z = norminv((1/2/m_z):(1/m_z):(1-1/2/m_z),0,sigma_z)';

disp('Solving the differential equations...');
tspan = linspace(0,12,10000);
init_cond = repmat([0.5,0.7,0,0,0,0],1,length(w));
if(~chance)
    r = norminv(0.9);
    xub = 100;
else
    r = norminv(0.95);
    xub = 6.6;
end
hl = @(x)(x(:,:,1)-xub);
dhldt0 = @(x)(x(:,:,1)-(1+theta(:,1))'.*x(:,:,1).*x(:,:,2));
dhldt1 = @(x)(0.4*x(:,:,1));
munf = @(x,n)((hl(x).^n)*w);
cf = @(x)(-2*(hl(x)*w).*(dhldt0(x)*w)+r^2/(1+r^2)*2*(hl(x).*dhldt0(x))*w)./...
    (-2*(hl(x)*w).*(dhldt1(x)*w)+r^2/(1+r^2)*2*(hl(x).*dhldt1(x))*w);
[t,x] = ode45(@(t,y)lveq(t,y,theta,Dtrain,munf,r,cf),tspan,init_cond,odeset('RelTol',1e-9,'AbsTol',1e-12));
x = reshape(x,length(tspan),6,[]);
for jj = 1:size(theta,1)
yyaxis right
if (jj<=4&&~betad&&fast)||(jj<=30&&~betad&&~fast)||(jj<=5&&betad)
l2=plot(t,x(:,2,jj),'-','LineWidth',3,'Color',custom_cols(1,:));
elseif (jj<=8&&jj>4&&~betad&&fast)||(jj<=60&&jj>30&&~betad&&~fast)||(jj<=11&&jj>5&&betad)
plot(t,x(:,2,jj),'-','LineWidth',1.75,'Color',custom_cols(1,:));
elseif (jj<=12&&jj>8&&~betad&&fast)||(jj<=90&&jj>60&&~betad&&~fast)||(jj<=19&&jj>11&&betad)
plot(t,x(:,2,jj),'-','LineWidth',0.5,'Color',custom_cols(1,:));
elseif (jj<=120&&jj>90&&~betad&&~fast)
plot(t,x(:,2,jj),':','LineWidth',0.5,'Color',custom_cols(1,:));
end
hold on
end
mu1 = munf(permute(x,[1,3,2]),1);
mu2 = munf(permute(x,[1,3,2]),2);
mu = -mu1.^2+r^2/(1+r^2)*mu2;
uopt = zeros(size(tspan));
for i = 1:length(tspan)
    if(mu(i)<=0)
        uopt(i) = arc_eval(tspan(i),Dtrain);
    else
        uopt(i) = max(0,cf(permute(x(i,:,:),[1,3,2])));
    end
end
if(chance)
yyaxis right
l3 = plot(t,mu1+xub,'-','LineWidth',3,'Color',custom_cols(2,:));
plot(t,mu1+xub+r*sqrt(max(0,mu2-mu1.^2)),'--',...
    t,mu1+xub-r*sqrt(max(0,mu2-mu1.^2)),'--','LineWidth',1.5,'Color',custom_cols(2,:));
end
yyaxis left
l1 = plot(t,uopt,'-','LineWidth',3,'Color',custom_cols(3,:));

%% Plot properties 
box 'on'
set(gca,'TickLabelInterpreter','latex');
ax = gca;
ax.FontSize = 15; 
ax.GridAlpha = 0.1;
ax.XAxis.LineWidth = 3;
ax.XAxis.Color = 'k';
ax.XAxis.Limits = [0,12];

ax.YAxis(1).LineWidth = 3;
ax.YAxis(1).Color= custom_cols(3,:);
ax.YAxis(1).Limits = [-0.01 1.01];

ax.YAxis(2).LineWidth = 3;
ax.YAxis(2).Color= custom_cols(1,:);
ax.YAxis(2).Limits = [-0.08 8.08];

hXL = xlabel('{$t$}','Interpreter','latex');
set(hXL,'FontSize',25);
set(ax,'xminorgrid','off','yminorgrid','off');
ph = [l1 l2];
LL = {'$u_1^{*}(t)$','$x_2^{*}(t;${\boldmath$\theta$}$)$'};
if(chance)
    ph = [ph,l3];
    LL = [LL,'$x_1^{*}(t;${\boldmath$\theta$}$)$'];
end
legend(ph,LL,'Interpreter','latex','FontSize',22,'LineWidth',1.5,'Location','northwest');

disp('Computing the expected utility and variances...');
u_d_1 = 0;
u_mc1_1 = 0;
u_mc2_1 = 0;
u_mc3_1 = 0;
for jj = 1:size(theta,1)
    u_mc1_2 = 0;
    u_mc2_2 = 0;
    u_mc3_2 = 0;
    for kk = 1:m_z
        u_mc2_3 = 0;
        u_mc3_3 = 0;
        for ll = 1:size(theta,1)
            u_mc2_3 = u_mc2_3+w(ll)*exp(-(Z(kk)+x(end,2,jj)-x(end,2,ll))^2/sigma_z^2/2)/((2*pi)^(1/2)*(sigma_z^2)^(1/2));
            u_mc3_3 = u_mc3_3+w(ll)*exp(-(((1+x(end,2,jj)/10)*Z(kk)+x(end,2,jj)-x(end,2,ll))/(1+x(end,2,ll)/10))^2/sigma_z^2/2)...
                /((2*pi)^(1/2)*(sigma_z^2)^(1/2))/(1+x(end,2,ll)/10)*(1+x(end,2,jj)/10);
        end
        u_mc1_2 = u_mc1_2+1/m_z*(-Z(kk)^2/sigma_z^2/2-log((2*pi)^(1/2)*(sigma_z^2)^(1/2))...
            +(Z(kk)+x(end,5:6,jj)*theta(jj,:)')^2/(sigma_z^2+x(end,5:6,jj)*K*x(end,5:6,jj)')/2+log((2*pi)^(1/2)*(sigma_z^2+x(end,5:6,jj)*K*x(end,5:6,jj)')^(1/2)));
        u_mc2_2 = u_mc2_2+1/m_z*(-Z(kk)^2/sigma_z^2/2-log((2*pi)^(1/2)*(sigma_z^2)^(1/2))-log(u_mc2_3));
        u_mc3_2 = u_mc3_2+1/m_z*(-Z(kk)^2/sigma_z^2/2-log((2*pi)^(1/2)*(sigma_z^2)^(1/2))-log(u_mc3_3));
    end
    u_d_1 = u_d_1+w(jj)*log(det(x(end,5:6,jj)'*sigma_z^(-2)*x(end,5:6,jj)+inv(K)));
    u_mc1_1 = u_mc1_1+w(jj)*u_mc1_2;
    u_mc2_1 = u_mc2_1+w(jj)*u_mc2_2;
    u_mc3_1 = u_mc3_1+w(jj)*u_mc3_2;
end
u_d = -u_d_1;
u_mc1 = -log(det(inv(K)))-2*u_mc1_1;
if(~betad)
u_mc2 = -log(det(inv(K)))-2*u_mc2_1;
u_mc3 = -log(det(inv(K)))-2*u_mc3_1;
else
u_mc2 = -u_mc2_1;
u_mc3 = -u_mc3_1;
end
gv = reshape(x(end,2,:),[],1);
Jv = reshape(1+x(end,2,:)/10,[],1);
m_pr = variance_y(0,1,theta,w,gv,Jv,sigma_z);
m_pi = integral(@(y)variance_y(y,2,theta,w,gv,Jv,sigma_z),-Inf,Inf,'ArrayValued',true,'AbsTol',1e-10,'RelTol',1e-6);
m_pd = integral(@(y)variance_y(y,3,theta,w,gv,Jv,sigma_z),-Inf,Inf,'ArrayValued',true,'AbsTol',1e-10,'RelTol',1e-6);
disp('Expected utility u_{D}:');
disp(u_d);
disp('Expected utility u_{MC} for state-independent Gaussian noise and normal prior:');
disp(u_mc1);
disp('Expected utility u_{MC} for state-independent noise:');
disp(u_mc2);
disp('Expected utility u_{MC} for state-dependent noise:');
disp(u_mc3);
disp('Covariance matrix for prior parameter distribution:');
disp(m_pr);
disp('Covariance matrix for posterior parameter distribution with state-independent noise:');
disp(m_pi);
disp('Covariance matrix for posterior parameter distribution with state-dependent noise:');
disp(m_pd);

return
