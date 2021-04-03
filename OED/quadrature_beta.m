%% Find points and weights for beta distributed priors 
%% Initialize
clear
close all
clc
profile on
nparam = 2;
% max_order = 7;
max_order = 9;
% max_order = 15;
a = [2;2];
b = [2;2];
c = [1/2;1/2];
orth_pols = @(k,n)sqrt((2*n+a(k)+b(k)-1)*gamma(n+a(k)+b(k)-1)*factorial(n)/gamma(n+a(k))/gamma(n+b(k))*beta(a(k),b(k)))*jacobiP(n,a(k)-1,b(k)-1,sym('x')/c(k));
pdf = @(x,k)betapdf((1-x/c(k))/2,a(k),b(k))/c(k)/2;

%% Choose number of training points 
% mtrain = ceil(nchoosek(max_order+nparam,nparam)/(nparam+1));
% mtrain = ceil(nchoosek(max_order+nparam,nparam));
% np = 10;
np = (max_order+1)/2;

%% Initial guess for theta
% rng(42);
% theta = c'.*(1-2*betarnd(ones(mtrain,1)*a',ones(mtrain,1)*b'));
% th1gv = c(1)*(1-2*betainv((1/2/np):(1/np):(1-1/2/np),a(1),b(1)));
% th2gv = c(2)*(1-2*betainv((1/2/np):(1/np):(1-1/2/np),a(2),b(2)));
% [th1,th2] = meshgrid(th1gv,th2gv);
% theta = [th1(:),th2(:)];
% theta = [mu(1),mu(2);mu(1)+2*sigma(1),mu(2);mu(1)-2*sigma(1),mu(2);...
%     mu(1)+sigma(1),mu(2)+sqrt(3)*sigma(2);mu(1)+sigma(1),mu(2)-sqrt(3)*sigma(2);...
%     mu(1)-sigma(1),mu(2)+sqrt(3)*sigma(2);mu(1)-sigma(1),mu(2)-sqrt(3)*sigma(2)];
% theta = [
%   -0.152631513018852  -0.152631513018852
%   -0.152631513018852   0.152631513018852
%    0.152631513018852  -0.152631513018852
%    0.152631513018852   0.152631513018852
%   -0.408248290463863   0.000000000000000
%    0.408248290463863   0.000000000000000
%    0.000000000000000  -0.408248290463863
%    0.000000000000000   0.408248290463863
%   -0.347631651864039  -0.347631651864039
%   -0.347631651864039   0.347631651864039
%    0.347631651864039  -0.347631651864039
%    0.347631651864039   0.347631651864039];
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
   0.440109933805545   0.210778895104562];
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
% w = [
%    0.161787835567326
%    0.161787835567326
%    0.161787835567326
%    0.161787835567326
%    0.051428571428571
%    0.051428571428571
%    0.051428571428571
%    0.051428571428571
%    0.036783593004103
%    0.036783593004103
%    0.036783593004103
%    0.036783593004103];
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
   0.016188049850202];
% w = [
%     0.0201051368113866*ones(30,1)
%     0.0119139564145941*ones(30,1)
%     0.0012962636171673*ones(30,1)
%     0.0000179764901854*ones(30,1)];
% wf = @(k,thetak)(2*np+a(k)+b(k)-1)./((c(k)^2-thetak.^2).*double(subs(diff(orth_pols(k,np),sym('x')),sym('x'),thetak)).^2);
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