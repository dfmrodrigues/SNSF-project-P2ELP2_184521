function [s,pr_inf] = solution_pyrrole_fA_fB(R,Ea,A,dH,cAin,cBin,V,T,fA,fB,noise)
k1 = A(1)*exp(-Ea(1)/R/T);
k2 = A(2)*exp(-Ea(2)/R/T);
g{1} = @(x)cAin*fA-x(1)*(fA+fB)/V-x(1)*x(2)*k1/V;
g{2} = @(x)cBin*fB-x(2)*(fA+fB)/V-x(1)*x(2)*k1/V-2*x(2)^2*k2/V;
dg{1,1} = @(x)-fA/V-fB/V-x(2)*k1/V;
dg{1,2} = @(x)-x(1)*k1/V;
dg{2,1} = @(x)-x(2)*k1/V;
dg{2,2} = @(x)-fA/V-fB/V-x(1)*k1/V-4*x(2)*k2/V;
s0 = [50;100];
s = newton_method(dg,g,s0,1000);
cA = s(1)/V;
cB = s(2)/V;
s(3) = V*k1*cA*cB*(cAin*(cBin-cB)-cA*cBin)/(cB*(cA*(cAin+cBin)*k1+2*cAin*cB*k2));
s(4) = V*k2*cB^2*(cAin*(cBin-cB)-cA*cBin)/(cB*(cA*(cAin+cBin)*k1+2*cAin*cB*k2));
s(5) = 0;
s(6) = -V*(-dH(1)*k1*cA*cB-dH(2)*k2*cB^2);
s(7) = cA*cB*((cA-cB+cBin)*k1+2*cB*k2)*V/(cAin*(cBin-cB)-cA*cBin);
s(8) = cB*(-cA^2*k1+cA*(cAin+cB)*k1-2*cA*cB*k2+2*cAin*cB*k2)*V/(cAin*(cBin-cB)-cA*cBin);
s(9) = V*k1*cA*cB;
s(10) = V*k2*cB^2;
s = normrnd(s,noise);
pr_inf = [];
end
