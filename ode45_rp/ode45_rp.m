function varargout = ode45_rp(ode,tspan,y0,rtol,atol,events,varargin)
%ODE45_RP  Solve non-stiff differential equations, medium order method.
%   [TOUT,YOUT] = ODE45_RP(ODE,TSPAN,Y0,RTOL,ATOL,EVENTS) integrates the
%   system of differential equations y' = f(t,y) from time T0 to TFINAL with
%   initial conditions Y0, and obtains solutions at specific times
%   TSPAN = [T0 T1 ... TFINAL] (all increasing or all decreasing). ODE is a
%   function handle, which results from the compilation of ode_mex.c, rp45.c
%   and ntrp45.c with C functions odefun and eventf (see definition later).
%   For a scalar T and a vector Y of size NEQ, odefun(F,T,Y,NEQ) must return
%   a column vector F corresponding to f(T,Y). Each column in the solution
%   array YOUT corresponds to a time returned in the row vector TOUT.
%
%   ODE45(ODEFUN,TSPAN,Y0,RTOL,ATOL,EVENTS) solves as above with the scalar
%   relative error tolerance RTOL and scalar absolute error tolerance ATOL.
%
%   [TOUT,YOUT,TE,YE,IE] = ODE45(ODE,TSPAN,Y0,RTOL,ATOL,EVENTS) with a
%   non-empty argument EVENTS, solves as above while also finding where
%   functions of (T,Y), called event functions, are zero, and terminates the
%   integration at a zero of one of these event functions. For each function
%   you specify whether the direction of the zero crossing matters. EVENTS
%   is a cell array with two elements: EVENTS{1} is a function handle that
%   results from the compilation of eventf_mex.c with a C function eventf,
%   and EVENTS{2} is the number of event functions NV. For a scalar T and a
%   vector Y of size NEQ, eventf(T,Y,V,DIRECTION,NV,NEQ) must return, for
%   the I-th event function: the value of the function V(I); DIRECTION(I)=0
%   if all zeros are to be computed (the default), +1 if only zeros where
%   the event function is increasing, and -1 if only zeros where the event
%   function is decreasing. Output TE is a row vector of times at which
%   events occur. Columns of YE are the corresponding solutions, and indices
%   in vector IE specify which event occurred.
%
%   Class support for inputs TSPAN, Y0
%     float: double

%   ODE45_RP is an implementation of the explicit Runge-Kutta (4,5) pair of
%   Dormand and Prince called variously RK5(4)7FM, DOPRI5, DP(4,5) and DP54.
%   It uses a "free" interpolant of order 4 communicated privately by
%   Dormand and Prince.  Local extrapolation is done.

%   Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
%   M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-14-94
%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 5.74.4.13 $  $Date: 2011/04/16 06:38:58 $

solver_name = 'ode45';

% Stats
nsteps  = 0;
nfailed = 0;
nfevals = 0;


% Handle solver arguments
if isempty(tspan) || isempty(y0)
    error(message('MATLAB:odearguments:TspanOrY0NotSupplied', solver_name));
end
if length(tspan) < 2
    error(message('MATLAB:odearguments:SizeTspan', solver_name));
end
htspan = abs(tspan(2) - tspan(1));
tspan = tspan(:);
ntspan = length(tspan);
t0 = tspan(1);
next = 2;       % next entry in tspan
tfinal = tspan(end);
args = varargin;                 % use f(t,y,p1,p2...)

y0 = y0(:);
neq = length(y0);

% Test that tspan is internally consistent.
if t0 == tfinal
    error(message('MATLAB:odearguments:TspanEndpointsNotDistinct'));
end
tdir = sign(tfinal - t0);
if any( tdir*diff(tspan) <= 0 )
    error(message('MATLAB:odearguments:TspanNotMonotonic'));
end

dataType = 'double';

if rtol < 100 * eps(dataType)
    rtol = 100 * eps(dataType);
    warning(message('MATLAB:odearguments:RelTolIncrease', sprintf( '%g', rtol )))
end
threshold = atol / rtol;

% By default, hmax is 1/10 of the interval.
hmax = abs(0.1*(tfinal-t0));
if hmax <= 0
    error(message('MATLAB:odearguments:MaxStepLEzero'));
end

% Handle the event function
teout = [];
yeout = [];
ieout = [];

if isempty(events)
    haveeventfun = 0;   % no Events function
    nv = 0;
    v0 = [];
else
    haveeventfun = 1;   % there is an Events function
    eventfun = events{1};
    nv = events{2};
    v0 = eventfun(t0,y0,neq,nv,args{:});
end

nfevals = nfevals + 1;

[tnew,ynew,h,f,err,nfevals,nfailed,nofailed,absh,nsteps,stop,...
    vnew,direction,indzc,v,nout,tout,yout,next,t,y] =...
    ode(t0,y0,neq,tfinal,htspan,threshold,nfevals,nfailed,rtol,hmax,tdir,nsteps,...
    v0,next,ntspan,haveeventfun,tspan,nv,args{:});

if haveeventfun && stop
    indzc = find(indzc);
    [te,ye,ie] = ...
        odezero(@ntrp45_mex,eventfun,args,v,t,y,vnew,tnew,ynew,neq,nv,t0,...
        tdir,direction,indzc,h,f);
    teout = [teout, te];
    yeout = [yeout, ye];
    ieout = [ieout, ie];
    
    tnew = te(end);
    ynew = ye(:,end);
    [oldnout,nout,tout_new,yout_new,next] =...
        rp45_mex(next,ntspan,haveeventfun,stop,tdir,tspan,tnew,ynew,t,y,h,f,...
        nofailed,err,rtol,absh,nout,neq);
    idx = oldnout+1:nout;
    tout(idx) = tout_new;
    yout(:,idx) = yout_new;
end

solver_output{1} = tout(1:nout);
solver_output{2} = yout(:,1:nout);
if haveeventfun
    solver_output{3} = teout;
    solver_output{4} = yeout;
    solver_output{5} = ieout;
end
varargout = solver_output;
end

function [te,ye,ie,stop] = ...
    odezero(ntrpfun,eventfun,eventargs,v,t,y,vnew,tnew,ynew,neq,nv,t0,...
    tdir,direction,indzc,h,f)
%ODEZERO Locate any zero-crossings of event functions in a time step.
%   ODEZERO is an event location helper function for the ODE Suite.  ODEZERO
%   uses Regula Falsi and information passed from the ODE solver to locate
%   any zeros in the half open time interval (T,TNEW] of the event functions
%   coded in eventfun.
%
%   See also ODE45, ODE23, ODE113, ODE15S, ODE23S, ODE23T, ODE23TB.

%   Mark W. Reichelt, Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.26.4.6 $  $Date: 2010/08/23 23:09:29 $

tol = 128*max(eps(t),eps(tnew));
tol = min(tol, abs(tnew - t));
te = [];
ye = [];
ie = [];
stop = 0;
rmin = realmin;

% Set up tL, tR, yL, yR, vL, vR, isterminal and direction.
tL = t;
yL = y;
vL = v;
tR = tnew;
yR = ynew;
vR = vnew;

% Initialize ttry so that we won't extrapolate if vL or vR is zero.
ttry = tR;

% Find all events before tnew or the first terminal event.
while true
    
    lastmoved = 0;
    while true
        % Events of interest shouldn't have disappeared, but new ones might
        % be found in other elements of the v vector.
        if isempty(indzc)
            if lastmoved ~= 0
                error(message('MATLAB:odezero:LostEvent'));
            end
            return;
        end
        
        % Check if the time interval is too short to continue looking.
        delta = tR - tL;
        if abs(delta) <= tol
            break;
        end
        
        if (tL == t) && any(vL(indzc) == 0 & vR(indzc) ~= 0)
            ttry = tL + tdir*0.5*tol;
            
        else
            % Compute Regula Falsi change, using leftmost possibility.
            change = 1;
            for j = indzc(:)'
                % If vL or vR is zero, try using old ttry to extrapolate.
                if vL(j) == 0
                    if (tdir*ttry > tdir*tR) && (vtry(j) ~= vR(j))
                        maybe = 1.0 - vR(j) * (ttry-tR) / ((vtry(j)-vR(j)) * delta);
                        if (maybe < 0) || (maybe > 1)
                            maybe = 0.5;
                        end
                    else
                        maybe = 0.5;
                    end
                elseif vR(j) == 0.0
                    if (tdir*ttry < tdir*tL) && (vtry(j) ~= vL(j))
                        maybe = vL(j) * (tL-ttry) / ((vtry(j)-vL(j)) * delta);
                        if (maybe < 0) || (maybe > 1)
                            maybe = 0.5;
                        end
                    else
                        maybe = 0.5;
                    end
                else
                    maybe = -vL(j) / (vR(j) - vL(j)); % Note vR(j) ~= vL(j).
                end
                if maybe < change
                    change = maybe;
                end
            end
            change = change * abs(delta);
            
            % Enforce minimum and maximum change.
            change = max(0.5*tol, min(change, abs(delta) - 0.5*tol));
            
            ttry = tL + tdir * change;
        end
        
        % Compute vtry.
        ytry = ntrpfun(ttry,t,y,h,f);
        vtry = eventfun(ttry,ytry,neq,nv,eventargs{:});
        
        % Check for any crossings between tL and ttry.
        indzc = find((sign(vL) ~= sign(vtry)) & (direction .* (vtry - vL) >= 0));
        if ~isempty(indzc)
            % Move right end of bracket leftward, remembering the old value.
            tswap = tR; tR = ttry; ttry = tswap;
            yswap = yR; yR = ytry; ytry = yswap;
            vswap = vR; vR = vtry; vtry = vswap;
            % Illinois method.  If we've moved leftward twice, halve
            % vL so we'll move closer next time.
            if lastmoved == 2
                % Watch out for underflow and signs disappearing.
                maybe = 0.5 * vL;
                i = find(abs(maybe) >= rmin);
                vL(i) = maybe(i);
            end
            lastmoved = 2;
        else
            % Move left end of bracket rightward, remembering the old value.
            tswap = tL; tL = ttry; ttry = tswap;
            yswap = yL; yL = ytry; ytry = yswap;
            vswap = vL; vL = vtry; vtry = vswap;
            % Illinois method.  If we've moved rightward twice, halve
            % vR so we'll move closer next time.
            if lastmoved == 1
                % Watch out for underflow and signs disappearing.
                maybe = 0.5 * vR;
                i = find(abs(maybe) >= rmin);
                vR(i) = maybe(i);
            end
            lastmoved = 1;
        end
        indzc = find((sign(vL) ~= sign(vR)) & (direction .* (vR - vL) >= 0));
    end
    
    j = ones(1,length(indzc));
    te = [te, tR(j)];
    ye = [ye, yR(:,j)];
    ie = [ie, indzc'];
    if any(indzc)
        if tL ~= t0
            stop = 1;
        end
        break;
    elseif abs(tnew - tR) <= tol
        %  We're not going to find events closer than tol.
        break;
    else
        % Shift bracket rightward from [tL tR] to [tR+0.5*tol tnew].
        ttry = tR; ytry = yR; vtry = vR;
        tL = tR + tdir*0.5*tol;
        yL = ntrpfun(tL,t,y,h,f);
        vL = eventfun(tL,yL,neq,nv,eventargs{:});
        tR = tnew; yR = ynew; vR = vnew;
    end
end
end
