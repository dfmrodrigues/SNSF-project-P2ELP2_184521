%% Find points and weights for normally distributed priors 
%% Initialize
clear
close all
clc
profile on
nparam = 2;
max_order = 7;
% max_order = 15;
mu  = [0;0];
sigma = [0.2;0.2];
orth_pols = @(k,n)hermiteH(n,(sym('x')-mu(k))/sqrt(2*sigma(k)^2))*(1/2)^(n/2)/sqrt(factorial(n));
pdf = @(x,k)normpdf(x,mu(k),sigma(k));

%% Choose number of training points 
% mtrain = ceil(nchoosek(max_order+nparam,nparam)/(nparam+1));
% mtrain = ceil(nchoosek(max_order+nparam,nparam));
% np = 10;
np = (max_order+1)/2;

%% Initial guess for theta
% rng(42);
% theta = normrnd(ones(mtrain,1)*mu',ones(mtrain,1)*sigma');
% th1gv = norminv((1/2/np):(1/np):(1-1/2/np),mu(1),sigma(1));
% th2gv = norminv((1/2/np):(1/np):(1-1/2/np),mu(2),sigma(2));
% [th1,th2] = meshgrid(th1gv,th2gv);
% theta = [th1(:),th2(:)];
% theta = [mu(1),mu(2);mu(1)+2*sigma(1),mu(2);mu(1)-2*sigma(1),mu(2);...
%     mu(1)+sigma(1),mu(2)+sqrt(3)*sigma(2);mu(1)+sigma(1),mu(2)-sqrt(3)*sigma(2);...
%     mu(1)-sigma(1),mu(2)+sqrt(3)*sigma(2);mu(1)-sigma(1),mu(2)-sqrt(3)*sigma(2)];
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
   0.396335765891742   0.396335765891742];
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
   0.210491910111102
   0.210491910111102
   0.210491910111102
   0.210491910111102
   0.027777777777778
   0.027777777777778
   0.027777777777778
   0.027777777777778
   0.011730312111120
   0.011730312111120
   0.011730312111120
   0.011730312111120];
% w = [
%     0.0201051368113866*ones(30,1)
%     0.0119139564145941*ones(30,1)
%     0.0012962636171673*ones(30,1)
%     0.0000179764901854*ones(30,1)];
% wf = @(k,thetak)(1/sigma(k)^2)./(double(subs(diff(orth_pols(k,np),sym('x')),sym('x'),thetak)).^2);
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