%% Find points and weights for gamma distributed priors 
%% Initialize
clear
close all
clc
profile on
nparam = 2;
max_order = 7;
% max_order = 15;
a = [2;2];
b = [1/25;1/25];
orth_pols = @(k,n)sqrt(factorial(n)*gamma(a(k))/gamma(n+a(k)))*laguerreL(n,a(k)-1,sym('x')/b(k));
pdf = @(x,k)gampdf(x,a(k),b(k));

%% Choose number of training points 
% mtrain = ceil(nchoosek(max_order+nparam,nparam)/(nparam+1));
% mtrain = ceil(nchoosek(max_order+nparam,nparam));
% np = 10;
np = (max_order+1)/2;

%% Initial guess for theta
% rng(42);
% theta = gamrnd(ones(mtrain,1)*a',ones(mtrain,1)*b');
% th1gv = gaminv((1/2/np):(1/np):(1-1/2/np),a(1),b(1));
% th2gv = gaminv((1/2/np):(1/np):(1-1/2/np),a(2),b(2));
% [th1,th2] = meshgrid(th1gv,th2gv);
% theta = [th1(:),th2(:)];
% theta = [mu(1),mu(2);mu(1)+2*sigma(1),mu(2);mu(1)-2*sigma(1),mu(2);...
%     mu(1)+sigma(1),mu(2)+sqrt(3)*sigma(2);mu(1)+sigma(1),mu(2)-sqrt(3)*sigma(2);...
%     mu(1)-sigma(1),mu(2)+sqrt(3)*sigma(2);mu(1)-sigma(1),mu(2)-sqrt(3)*sigma(2)];
theta = [
   0.227893614292822   0.035807733194508
   0.441601349350130   0.199971780204322
   0.104646413230098   0.083444284794823
   0.030672982142861   0.118781141172585
   0.026345261608404   0.266670653710175
   0.437554555452437   0.052183711420736
   0.068350690046067   0.462563746562234
   0.100072891017704   0.020046564744642
   0.231916162975243   0.134238647128262
   0.029190643801748   0.034576159548474
   0.216396090616353   0.341393562675628
   0.101399066135641   0.212495067648374];
% theta = [
%     0.160635659707143*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]';...
%     0.373712306584429*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]';...
%     0.602436406397929*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]';...
%     0.866951943872305*[cos((2*pi/30):(2*pi/30):(2*pi));sin((2*pi/30):(2*pi/30):(2*pi))]'];
% thf = @(k)roots(double(coeffs(orth_pols(k,np),sym('x'),'All')));
% th1gv = thf(1);
% th2gv = thf(2);
% [th1,th2] = meshgrid(th1gv,th2gv);
% theta = [th1(:),th2(:)];

%% Initial guess for weights
% w = ones(mtrain,1)/mtrain;
% w = ones(np^2,1)/(np^2);
w = [
   0.043963291448661
   0.000230436936686
   0.277862738201366
   0.197376444960076
   0.016058351811116
   0.001081207941655
   0.000704023916693
   0.143747917037863
   0.028844376047048
   0.233164295534674
   0.001271321141883
   0.055695595022278];
% w = [
%     0.0201051368113866*ones(30,1)
%     0.0119139564145941*ones(30,1)
%     0.0012962636171673*ones(30,1)
%     0.0000179764901854*ones(30,1)];
% wf = @(k,thetak)(1/b(k))./(thetak.*double(subs(diff(orth_pols(k,np),sym('x')),sym('x'),thetak)).^2);
% w1gv = wf(1,th1gv);
% w2gv = wf(2,th2gv);
% [w1,w2] = meshgrid(w1gv,w2gv);
% w = w1(:).*w2(:);

%% Execute general code to find points and weights for any prior
[wopt,thetaopt,pols,dpols,innerp,diff_pols_eval,diff_pols_deval,diff_basis_eval,diff_basis_deval] = quadrature(nparam,max_order,orth_pols,pdf,w,theta);
disp(sum(sum(abs(w-wopt))))
disp(sum(sum(abs(theta-thetaopt))))
disp(sum(sum(sum(abs(innerp-eye(max_order+1))))))
disp(sum(sum(sum(abs(diff_pols_eval)))))
disp(sum(sum(sum(abs(diff_pols_deval)))))
disp(sum(sum(sum(abs(diff_basis_eval)))))
disp(sum(sum(sum(abs(diff_basis_deval)))))
profile viewer