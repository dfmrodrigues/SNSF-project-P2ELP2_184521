function dy_dt = dynamics_pyrrole_1(t,y,Z,qex,fA,fB,cAin,cBin,V,k10,k20,dH1,dH2,Ea1,Ea2,R,nCp,Tref,C,S,ck1,bk1,lags,eu,ey,h,u0)
Q = y(1);
n = y(2:5);
nA = n(1);
nB = n(2);
nC = n(3);
nD = n(4);
iqex = y(6);
ifA = y(7);
ifB = y(8);
T = Tref+Q/nCp;
r(1,:) = k10*exp(Ea1/R*(1/Tref-1/T))*nA*nB/V;
r(2,:) = k20*exp(Ea2/R*(1/Tref-1/T))*nB^2/V;
prevQ = [Z(1,:),Q];
prevn = [Z(2:5,:),n];
prevnA = prevn(1,:);
prevnB = prevn(2,:);
previqex = [Z(6,:),iqex];
previfA = [Z(7,:),ifA];
previfB = [Z(8,:),ifB];
taui = (t-[lags(1),lags])/h;
tauf = (t-[lags(2:end),0,0])/h;
dtau = tauf-taui;
flags = max(0,min(1,tauf./dtau));
prevqex = u0(1)*(1-flags)+(flags>0).*diff([previqex(1),previqex(1:end-1);previqex(2:end),previqex(end)])./(h*dtau);
prevfA = u0(2)*(1-flags)+(flags>0).*diff([previfA(1),previfA(1:end-1);previfA(2:end),previfA(end)])./(h*dtau);
prevfB = u0(3)*(1-flags)+(flags>0).*diff([previfB(1),previfB(1:end-1);previfB(2:end),previfB(end)])./(h*dtau);
dprevy = ey(t);
dprevQ = dprevy(1,:);
dprevn = dprevy(2:5,:);
dprevnA = dprevn(1,:);
dprevnB = dprevn(2,:);
dprevu = eu(t);
dprevqex = dprevu(1,:);
dprevfA = dprevu(2,:);
dprevfB = dprevu(3,:);
prevsa(1,:) = (prevqex-dprevqex)-((prevfA-dprevfA)+(prevfB-dprevfB)).*(prevQ+dprevQ)/V;
prevsa(2,:) = (prevfA-dprevfA)*cAin-((prevfA-dprevfA)+(prevfB-dprevfB)).*(prevnA+dprevnA)/V;
prevsa(3,:) = (prevfB-dprevfB)*cBin-((prevfA-dprevfA)+(prevfB-dprevfB)).*(prevnB+dprevnB)/V;
hatyu1 = S*C*[prevQ+dprevQ;prevn+dprevn]*ck1'-prevsa*bk1';
dQ = dprevQ(:,end);
dn = dprevn(:,end);
dqex = dprevqex(:,end);
dfA = dprevfA(:,end);
dfB = dprevfB(:,end);
qex = qex(t,Q+dQ,n+dn,hatyu1)+dqex;
fA = fA(t,Q+dQ,n+dn,hatyu1)+dfA;
fB = fB(t,Q+dQ,n+dn,hatyu1)+dfB;
dy_dt = [-dH1*r(1)-dH2*r(2)+qex-(fA+fB)*Q/V;...
    -r(1)+fA*cAin-(fA+fB)*nA/V;...
    -r(1)-2*r(2)+fB*cBin-(fA+fB)*nB/V;...
    r(1)-(fA+fB)*nC/V;...
    r(2)-(fA+fB)*nD/V;...
    qex;...
    fA;...
    fB];
end