function dy_dt = dynamics_appj(t,y,Z,P,Ba,Tinf,Tb,Tmax,Rho,Cp,d,r,eta,K,m1,m2,C,S,ck1,bk1,lags,eu,ey,h,T0)
% Gather current and past information
T = y(2); % Current temperature
flagT = y(3); % Current flag for switching to constraint-seeking arc
itP = y(4); % Current integral of applied power \tilde{P}
idT = y(5); % Current integral of temperature noise
prevT = [Z(2,:),T]; % Temperature in the filter window
previtP = [Z(4,:),itP]; % Integral of applied power \tilde{P} in the filter window
previdT = [Z(5,:),idT]; % Integral of temperature noise in the filter window
% Define constants for the plant
k1 = m2*eta/Rho/Cp/pi/r^2/d;
k2 = 2*pi*r*d*K*m1/Rho/Cp/pi/r^2/d*(1/2+(Tmax-Tb)/(Tb-Tinf))*log(1+(Tb-Tinf)/(Tmax-Tb));
P0 = k2*((T0-Tinf)-(T0-Tb))/log((T0-Tinf)/(T0-Tb))/k1; % Initial value of the manipulated variable P
% Implement FIR filter and unknown rate estimation
taui = (t-[lags(1),lags])/h;
tauf = (t-[lags(2:end),0,0])/h;
dtau = tauf-taui;
flags = max(0,min(1,tauf./dtau));
prevtP = P0*(1-flags)+(flags>0).*diff([previtP(1),previtP(1:end-1);previtP(2:end),previtP(end)])./(h*dtau); % Applied power \tilde{P} in the filter window
prevdT = 0*(1-flags)+(flags>0).*diff([previdT(1),previdT(1:end-1);previdT(2:end),previdT(end)])./(h*dtau); % Temperature noise in the filter window
prevsa(1,:) = Ba(prevT+prevdT).*prevtP; % Estimated known rate y_a
hatyu1 = S*C*(prevT+prevdT)*ck1'-prevsa*bk1'; % Estimated unknown rate y_u
% Gather current input and output disturbances
du = eu(t);
dP = du(1); % Input disturbance
dy = ey(t);
dT = dy(1); % Output disturbance
tP = P(t,T+dT,flagT,hatyu1); % Applied power \tilde{P} according to the specific control law
P = tP+dP; % True power P
% Conpute state dynamics (CEM, temperature, flag for switching to constraint-seeking arc, integral of applied power \tilde{P}, integral of temperature noise)
dy_dt = [0.5^(43.0+273.15-T)/60.0;...
    k1*P-k2*((T-Tinf)-(T-Tb))/log((T-Tinf)/(T-Tb));...
    1+erf((T-Tmax)/0.05);...
    tP;...
    dT];
end