function [tau,p_sol,d_sol,deg] = hessian_optim(n,N,p_mon,p_s,d,cphi,u0,sc)
ab_mon = p_mon;
a = [1;zeros(p_s-1,1)];
b = cphi{1};
c = round(-cell2mat(cphi(2:end)),10);
c = c(:,sum(abs(c)>1e-8,1)~=0);
w = 1;
idx = 1;
rad = sum((0.5*ones(1,N)).^n)^(1/n);
deg = n/2+1;
Ik = 1:size(ab_mon,2);
Kj = ones(1,size(c,2)+size(d,2));
deg_max = 10;
[tau,p_sol,d_sol,deg] = multiv_pols_min(ab_mon,a,b,[c,d],w,idx,rad,deg,Ik,Kj,deg_max);
% deg = n/2+1;
% Kj = ones(1,size(c,2));
% [tau,p_sol,d_sol,deg,si] = multiv_pols_min(ab_mon,a,b,c,w,idx,rad,deg,Ik,Kj,deg_max);
% objPoly.typeCone = 1;
% objPoly.dimVar = N;
% objPoly.degree = n;
% objPoly.noTerms = p_s;
% objPoly.supports = p_mon;
% objPoly.coef = b;
% ineqPolySys = cell(1,size(c,2)+1);
% for i = 1:size(c,2)
%     ineqPolySys{i}.typeCone = 1;
%     ineqPolySys{i}.dimVar = N;
%     ineqPolySys{i}.degree = n;
%     ineqPolySys{i}.noTerms = p_s;
%     ineqPolySys{i}.supports = p_mon;
%     ineqPolySys{i}.coef = c(:,i);
% end
% ineqPolySys{end}.typeCone = 1;
% ineqPolySys{end}.dimVar = N;
% ineqPolySys{end}.degree = n;
% ineqPolySys{end}.noTerms = N+1;
% ineqPolySys{end}.supports = [zeros(1,N);n*eye(N)];
% ineqPolySys{end}.coef = [rad^n;-ones(N,1)];
% lbd = -1e10*ones(1,N);
% ubd = 1e10*ones(1,N);
% param.relaxOrder = n/2+1;
% param.SDPsolverOutFile = 1;
% param.SDPsolver = 'sedumi';
% [param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);
% POP.objValue
% POP.xVect'
% param.relaxOrder
% tau-POP.objValue
% p_sol-POP.xVect'
% d_sol-POP.xVect'
% deg-param.relaxOrder
% si.SDP.A+SDPinfo.SeDuMiA
% si.SDP.b+SDPinfo.SeDuMib
% si.SDP.c-SDPinfo.SeDuMic
% 
% mset clear
% mpol x 5
% g0 = prod(((ones(size(b))*x').^ab_mon)')*b;
% c_mon = calc_mon(n,N);
% gi = [(prod(((ones(size(c,1),1)*x').^c_mon)')*c)' >= 0;...
%     sum(x.^n) <= rad^n];
% deg = n/2+1;
% status = 0;
% while(deg<100)
%     P = msdp(min(g0),gi,deg);
%     [status,obj,M] = msol(P);
%     if(status==1)
%         break;
%     end
%     deg = deg+1;
% end
% obj
% sol = squeeze(double(M))'
% deg
p_sol = u0+p_sol.*sc;
d_sol = u0+d_sol.*sc;
end