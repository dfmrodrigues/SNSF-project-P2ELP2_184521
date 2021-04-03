function varargout = ode45f(ode,tspan,y0,options,varargin)
%ODE45  Solve non-stiff differential equations, medium order method.
%   [TOUT,YOUT] = ODE45(ODEFUN,TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates
%   the system of differential equations y' = f(t,y) from time T0 to TFINAL 
%   with initial conditions Y0. ODEFUN is a function handle. For a scalar T
%   and a vector Y, ODEFUN(T,Y) must return a column vector corresponding 
%   to f(t,y). Each row in the solution array YOUT corresponds to a time 
%   returned in the column vector TOUT.  To obtain solutions at specific 
%   times T0,T1,...,TFINAL (all increasing or all decreasing), use TSPAN = 
%   [T0 T1 ... TFINAL].     
%   
%   [TOUT,YOUT] = ODE45(ODEFUN,TSPAN,Y0,OPTIONS) solves as above with default
%   integration properties replaced by values in OPTIONS, an argument created
%   with the ODESET function. See ODESET for details. Commonly used options 
%   are scalar relative error tolerance 'RelTol' (1e-3 by default) and vector
%   of absolute error tolerances 'AbsTol' (all components 1e-6 by default).
%   If certain components of the solution must be non-negative, use
%   ODESET to set the 'NonNegative' property to the indices of these
%   components.
%   
%   ODE45 can solve problems M(t,y)*y' = f(t,y) with mass matrix M that is
%   nonsingular. Use ODESET to set the 'Mass' property to a function handle 
%   MASS if MASS(T,Y) returns the value of the mass matrix. If the mass matrix 
%   is constant, the matrix can be used as the value of the 'Mass' option. If
%   the mass matrix does not depend on the state variable Y and the function
%   MASS is to be called with one input argument T, set 'MStateDependence' to
%   'none'. ODE15S and ODE23T can solve problems with singular mass matrices.  
%
%   [TOUT,YOUT,TE,YE,IE] = ODE45(ODEFUN,TSPAN,Y0,OPTIONS) with the 'Events'
%   property in OPTIONS set to a function handle EVENTS, solves as above 
%   while also finding where functions of (T,Y), called event functions, 
%   are zero. For each function you specify whether the integration is 
%   to terminate at a zero and whether the direction of the zero crossing 
%   matters. These are the three column vectors returned by EVENTS: 
%   [VALUE,ISTERMINAL,DIRECTION] = EVENTS(T,Y). For the I-th event function: 
%   VALUE(I) is the value of the function, ISTERMINAL(I)=1 if the integration 
%   is to terminate at a zero of this event function and 0 otherwise. 
%   DIRECTION(I)=0 if all zeros are to be computed (the default), +1 if only 
%   zeros where the event function is increasing, and -1 if only zeros where 
%   the event function is decreasing. Output TE is a column vector of times 
%   at which events occur. Rows of YE are the corresponding solutions, and 
%   indices in vector IE specify which event occurred.    
%
%   SOL = ODE45(ODEFUN,[T0 TFINAL],Y0...) returns a structure that can be
%   used with DEVAL to evaluate the solution or its first derivative at 
%   any point between T0 and TFINAL. The steps chosen by ODE45 are returned 
%   in a row vector SOL.x.  For each I, the column SOL.y(:,I) contains 
%   the solution at SOL.x(I). If events were detected, SOL.xe is a row vector 
%   of points at which events occurred. Columns of SOL.ye are the corresponding 
%   solutions, and indices in vector SOL.ie specify which event occurred. 
%
%   Example    
%         [t,y]=ode45(@vdp1,[0 20],[2 0]);   
%         plot(t,y(:,1));
%     solves the system y' = vdp1(t,y), using the default relative error
%     tolerance 1e-3 and the default absolute tolerance of 1e-6 for each
%     component, and plots the first component of the solution. 
%   
%   Class support for inputs TSPAN, Y0, and the result of ODEFUN(T,Y):
%     float: double, single
%
%   See also ODE23, ODE113, ODE15S, ODE23S, ODE23T, ODE23TB, ODE15I,
%            ODESET, ODEPLOT, ODEPHAS2, ODEPHAS3, ODEPRINT, DEVAL,
%            ODEEXAMPLES, RIGIDODE, BALLODE, ORBITODE, FUNCTION_HANDLE.

%   ODE45 is an implementation of the explicit Runge-Kutta (4,5) pair of
%   Dormand and Prince called variously RK5(4)7FM, DOPRI5, DP(4,5) and DP54.
%   It uses a "free" interpolant of order 4 communicated privately by
%   Dormand and Prince.  Local extrapolation is done.

%   Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
%   M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-14-94
%   Copyright 1984-2017 The MathWorks, Inc.

solver_name = 'ode45';

% Check inputs
if nargin < 4
  options = [];
  if nargin < 3
    y0 = [];
    if nargin < 2
      tspan = [];
      if nargin < 1
        error(message('MATLAB:ode45:NotEnoughInputs'));
      end  
    end
  end
end

% Stats
nsteps  = 0;
nfailed = 0;
nfevals = 0; 

% Output
FcnHandlesUsed  = isa(ode,'function_handle');
output_sol = (FcnHandlesUsed && (nargout==1));      % sol = odeXX(...)
output_ty  = (~output_sol && (nargout > 0));  % [t,y,...] = odeXX(...)
% There might be no output requested...

sol = []; f3d = [];
if output_sol
  sol.solver = solver_name;
  sol.extdata.odefun = ode;
  sol.extdata.options = options;
  sol.extdata.varargin = varargin;
end

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, odeFcn, ...
  options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, dataType] = ...
  odearguments(FcnHandlesUsed, solver_name, ode, tspan, y0, options, varargin);
nfevals = nfevals + 1;

% Handle the output
if nargout > 0
  outputFcn = odeget(options,'OutputFcn',[],'fast');
else
  outputFcn = odeget(options,'OutputFcn',@odeplot,'fast');
end
outputArgs = {};
if isempty(outputFcn)
  haveOutputFcn = false;
else
  haveOutputFcn = true;
  outputs = odeget(options,'OutputSel',1:neq,'fast');
  if isa(outputFcn,'function_handle')
    % With MATLAB 6 syntax pass additional input arguments to outputFcn.
    outputArgs = varargin;
  end
end
refine = max(1,odeget(options,'Refine',4,'fast'));
if ntspan > 2
  outputAt = 1;          % output only at tspan points
elseif refine <= 1
  outputAt = 2;          % computed points, no refinement
else
  outputAt = 3;          % computed points, with refinement
  S = (1:refine-1) / refine;
end
printstats = strcmp(odeget(options,'Stats','off','fast'),'on');

% Handle the event function
[haveEventFcn,eventFcn,eventArgs,valt,teout,yeout,ieout] = ...
  odeevents(FcnHandlesUsed,odeFcn,t0,y0,options,varargin);

% Handle the mass matrix
[Mtype, M, Mfun] =  odemass(FcnHandlesUsed,odeFcn,t0,y0,options,varargin);
if Mtype > 0  % non-trivial mass matrix
  Msingular = odeget(options,'MassSingular','no','fast');
  if strcmp(Msingular,'maybe')
    warning(message('MATLAB:ode45:MassSingularAssumedNo'));
  elseif strcmp(Msingular,'yes')
    error(message('MATLAB:ode45:MassSingularYes'));
  end
  % Incorporate the mass matrix into odeFcn and odeArgs.
  [odeFcn,odeArgs] = odemassexplicit(FcnHandlesUsed,Mtype,odeFcn,odeArgs,Mfun,M);
  f0 = feval(odeFcn,t0,y0,odeArgs{:});
  nfevals = nfevals + 1;
end

% Non-negative solution components
idxNonNegative = odeget(options,'NonNegative',[],'fast');
nonNegative = ~isempty(idxNonNegative);
if nonNegative  % modify the derivative function
  [odeFcn,thresholdNonNegative] = odenonnegative(odeFcn,y0,threshold,idxNonNegative);
  f0 = feval(odeFcn,t0,y0,odeArgs{:});
  nfevals = nfevals + 1;
end

t = t0;
y = y0;

% Allocate memory if we're generating output.
nout = 0;
tout = []; yout = [];
if nargout > 0
  if output_sol
    chunk = min(max(100,50*refine), refine+floor((2^11)/neq));
    tout = zeros(1,chunk,dataType);
    yout = zeros(neq,chunk,dataType);
    f3d  = zeros(neq,7,chunk,dataType);
  else
    if ntspan > 2                         % output only at tspan points
      tout = zeros(1,ntspan,dataType);
      yout = zeros(neq,ntspan,dataType);
    else                                  % alloc in chunks
      chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
      tout = zeros(1,chunk,dataType);
      yout = zeros(neq,chunk,dataType);
    end
  end
  nout = 1;
  tout(nout) = t;
  yout(:,nout) = y;
end

% Initialize method parameters.
pow = 1/5;
A = [1/5, 3/10, 4/5, 8/9, 1, 1]; % Still used by restarting criteria
% B = [
%     1/5         3/40    44/45   19372/6561      9017/3168       35/384
%     0           9/40    -56/15  -25360/2187     -355/33         0
%     0           0       32/9    64448/6561      46732/5247      500/1113
%     0           0       0       -212/729        49/176          125/192
%     0           0       0       0               -5103/18656     -2187/6784
%     0           0       0       0               0               11/84
%     0           0       0       0               0               0
%     ];
% E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];

% Same values as above extracted as scalars (1 and 0 values ommitted)
a2=cast(1/5,dataType);
a3=cast(3/10,dataType);
a4=cast(4/5,dataType);
a5=cast(8/9,dataType);

b11=cast(1/5,dataType); 
b21=cast(3/40,dataType); 
b31=cast(44/45,dataType);
b41=cast(19372/6561,dataType);
b51=cast(9017/3168,dataType);
b61=cast(35/384,dataType);
b22=cast(9/40,dataType);
b32=cast(-56/15,dataType);
b42=cast(-25360/2187,dataType);
b52=cast(-355/33,dataType);
b33=cast(32/9,dataType);
b43=cast(64448/6561,dataType);
b53=cast(46732/5247,dataType);
b63=cast(500/1113,dataType);
b44=cast(-212/729,dataType);
b54=cast(49/176,dataType);
b64=cast(125/192,dataType);
b55=cast(-5103/18656,dataType);
b65=cast(-2187/6784,dataType);
b66=cast(11/84,dataType);

e1=cast(71/57600,dataType);
e3=cast(-71/16695,dataType);
e4=cast(71/1920,dataType);
e5=cast(-17253/339200,dataType);
e6=cast(22/525,dataType);
e7=cast(-1/40,dataType);

hmin = 16*eps(t);
if isempty(htry)
  % Compute an initial step size h using y'(t).
  absh = min(hmax, htspan);
  if normcontrol
    rh = (norm(f0) / max(normy,threshold)) / (0.8 * rtol^pow);
  else
    rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^pow);
  end
  if absh * rh > 1
    absh = 1 / rh;
  end
  absh = max(absh, hmin);
else
  absh = min(hmax, max(hmin, htry));
end
f1 = f0;

% Initialize the output function.
if haveOutputFcn
  feval(outputFcn,[t tfinal],y(outputs),'init',outputArgs{:});
end

% Cleanup the main ode function call 
FcnUsed = isa(odeFcn,'function_handle');
odeFcn_main = odefcncleanup(FcnUsed,odeFcn,odeArgs);

% THE MAIN LOOP

done = false;
while ~done
  
  % By default, hmin is a small number such that t+hmin is only slightly
  % different than t.  It might be 0 if t is 0.
  hmin = 16*eps(t);
  absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
  h = tdir * absh;
  
  % Stretch the step if within 10% of tfinal-t.
  if 1.1*absh >= abs(tfinal - t)
    h = tfinal - t;
    absh = abs(h);
    done = true;
  end
  
  % LOOP FOR ADVANCING ONE STEP.
  nofailed = true;                      % no failed attempts
  while true
    y2 = y + h .* (b11.*f1 );
    t2 = t + h .* a2;
    f2 = odeFcn_main(t2, y2);
        
    y3 = y + h .* (b21.*f1 + b22.*f2 );
    t3 = t + h .* a3;
    f3 = odeFcn_main(t3, y3);
        
    y4 = y + h .* (b31.*f1 + b32.*f2 + b33.*f3 );
    t4 = t + h .* a4;
    f4 = odeFcn_main(t4, y4);
        
    y5 = y + h .* (b41.*f1 + b42.*f2 + b43.*f3 + b44.*f4 );
    t5 = t + h .* a5;
    f5 = odeFcn_main(t5, y5);
       
    y6 = y + h .* (b51.*f1 + b52.*f2 + b53.*f3 + b54.*f4 + b55.*f5 );
    t6 = t + h;
    f6 = odeFcn_main(t6, y6);

    tnew = t + h;
    if done
      tnew = tfinal;   % Hit end point exactly.
    end
    h = tnew - t;      % Purify h.     
    
    ynew = y + h.* ( b61.*f1 + b63.*f3 + b64.*f4 + b65.*f5 + b66.*f6 );
    f7 = odeFcn_main(tnew,ynew);
    
    nfevals = nfevals + 6;              
    
    % Estimate the error.
    NNrejectStep = false;
    fE = f1*e1 + f3*e3 + f4*e4 + f5*e5 + f6*e6 + f7*e7;
    if normcontrol
      normynew = norm(ynew);
      errwt = max(max(normy,normynew),threshold);
      err = absh * (norm(fE) / errwt);
      if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
        errNN = norm( max(0,-ynew(idxNonNegative)) ) / errwt ;
        if errNN > rtol
          err = errNN;
          NNrejectStep = true;
        end
      end      
    else
      err = absh * norm((fE) ./ max(max(abs(y),abs(ynew)),threshold),inf);      
      if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
        errNN = norm( max(0,-ynew(idxNonNegative)) ./ thresholdNonNegative, inf);      
        if errNN > rtol
          err = errNN;
          NNrejectStep = true;
        end
      end            
    end
    
    % Accept the solution only if the weighted error is no more than the
    % tolerance rtol.  Estimate an h that will yield an error of rtol on
    % the next step or the next try at taking this step, as the case may be,
    % and use 0.8 of this value to avoid failures.
    if err > rtol                       % Failed step
      nfailed = nfailed + 1;            
      if absh <= hmin
        warning(message('MATLAB:ode45:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));
        solver_output = odefinalize(solver_name, sol,...
                                    outputFcn, outputArgs,...
                                    printstats, [nsteps, nfailed, nfevals],...
                                    nout, tout, yout,...
                                    haveEventFcn, teout, yeout, ieout,...
                                    {f3d,idxNonNegative});
        if nargout > 0
          varargout = solver_output;
        end  
        return;
      end
      
      if nofailed
        nofailed = false;
        if NNrejectStep
          absh = max(hmin, 0.5*absh);
        else
          absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
        end
      else
        absh = max(hmin, 0.5 * absh);
      end
      h = tdir * absh;
      done = false;
      
    else                                % Successful step

      NNreset_f7 = false;
      if nonNegative && any(ynew(idxNonNegative)<0)
        ynew(idxNonNegative) = max(ynew(idxNonNegative),0);
        if normcontrol
          normynew = norm(ynew);
        end
        NNreset_f7 = true;
      end  
                  
      break;
      
    end
  end
  nsteps = nsteps + 1;                  
  
  if haveEventFcn
    f = [f1 f2 f3 f4 f5 f6 f7];
    [te,ye,ie,valt,stop] = ...
        odezero(@ntrp45,eventFcn,eventArgs,valt,t,y,tnew,ynew,t0,h,f,idxNonNegative);
    if ~isempty(te)
      if output_sol || (nargout > 2)
        teout = [teout, te];
        yeout = [yeout, ye];
        ieout = [ieout, ie];
      end
      if stop               % Stop on a terminal event.               
        % Adjust the interpolation data to [t te(end)].   
        
        % Update the derivatives using the interpolating polynomial.
        taux = t + (te(end) - t)*A;        
        [~,f(:,2:7)] = ntrp45(taux,t,y,[],[],h,f,idxNonNegative);
        f2 = f(:,2); f3 = f(:,3); f4 = f(:,4); f5 = f(:,5); f6 = f(:,6); f7 = f(:,7);
        
        tnew = te(end);
        ynew = ye(:,end);
        h = tnew - t;
        done = true;
      end
    end
  end

  if output_sol
    nout = nout + 1;
    if nout > length(tout)
      tout = [tout, zeros(1,chunk,dataType)];  % requires chunk >= refine
      yout = [yout, zeros(neq,chunk,dataType)];
      f3d  = cat(3,f3d,zeros(neq,7,chunk,dataType)); 
    end
    tout(nout) = tnew;
    yout(:,nout) = ynew;
    f3d(:,:,nout) = [f1 f2 f3 f4 f5 f6 f7];
  end  
    
  if output_ty || haveOutputFcn 
    switch outputAt
     case 2      % computed points, no refinement
      nout_new = 1;
      tout_new = tnew;
      yout_new = ynew;
     case 3      % computed points, with refinement
      tref = t + (tnew-t)*S;
      nout_new = refine;
      tout_new = [tref, tnew];
      yntrp45 = ntrp45split(tref,t,y,h,f1,f3,f4,f5,f6,f7,idxNonNegative);
      yout_new = [yntrp45, ynew];
     case 1      % output only at tspan points
      [oldnout,nout,tout_new,yout_new,next] = rp45_mex(next,ntspan,haveEventFcn,0,tdir,tspan,tnew,ynew,t,y,h,[f1,f2,f3,f4,f5,f6,f7],nofailed,err,rtol,absh,nout,neq);
      nout_new = nout-oldnout;
%       nout_new =  0;
%       tout_new = [];
%       yout_new = [];
%       while next <= ntspan  
%         if tdir * (tnew - tspan(next)) < 0
%           if haveEventFcn && stop     % output tstop,ystop
%             nout_new = nout_new + 1;
%             tout_new = [tout_new, tnew];
%             yout_new = [yout_new, ynew];            
%           end
%           break;
%         end
%         nout_new = nout_new + 1;              
%         tout_new = [tout_new, tspan(next)];
%         if tspan(next) == tnew
%           yout_new = [yout_new, ynew];            
%         else
%           yntrp45 = ntrp45split(tspan(next),t,y,h,f1,f3,f4,f5,f6,f7,idxNonNegative);
%           yout_new = [yout_new, yntrp45];
%         end  
%         next = next + 1;
%       end
    end
    
    if nout_new > 0
      if output_ty
%         oldnout = nout;
%         nout = nout + nout_new;
%         if nout > length(tout)
%           tout = [tout, zeros(1,chunk,dataType)];  % requires chunk >= refine
%           yout = [yout, zeros(neq,chunk,dataType)];
%         end
        idx = oldnout+1:nout;        
        tout(idx) = tout_new;
        yout(:,idx) = yout_new;
      end
      if haveOutputFcn
        stop = feval(outputFcn,tout_new,yout_new(outputs,:),'',outputArgs{:});
        if stop
          done = true;
        end  
      end     
    end  
  end
  
  if done
    break
  end

  % If there were no failures compute a new h.
  if nofailed
    % Note that absh may shrink by 0.8, and that err may be 0.
    temp = 1.25*(err/rtol)^pow;
    if temp > 0.2
      absh = absh / temp;
    else
      absh = 5.0*absh;
    end
  end
  
  % Advance the integration one step.
  t = tnew;
  y = ynew;
  if normcontrol
    normy = normynew;
  end
  if NNreset_f7
    % Used f7 for unperturbed solution to interpolate.  
    % Now reset f7 to move along constraint. 
    f7 = odeFcn_main(tnew,ynew);
    nfevals = nfevals + 1;
  end
  f1 = f7;  % Already have f(tnew,ynew)
  
end

solver_output = odefinalize(solver_name, sol,...
                            outputFcn, outputArgs,...
                            printstats, [nsteps, nfailed, nfevals],...
                            nout, tout, yout,...
                            haveEventFcn, teout, yeout, ieout,...
                            {f3d,idxNonNegative});
if nargout > 0
  varargout = solver_output;
end  
end

function [neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, args, odeFcn, ...
          options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, ...
          dataType ] =   ...
    odearguments(FcnHandlesUsed, solver, ode, tspan, y0, options, extras)
%ODEARGUMENTS  Helper function that processes arguments for all ODE solvers.
%
%   See also ODE113, ODE15I, ODE15S, ODE23, ODE23S, ODE23T, ODE23TB, ODE45.

%   Mike Karr, Jacek Kierzenka
%   Copyright 1984-2017 The MathWorks, Inc.

if strcmp(solver,'ode15i')
  FcnHandlesUsed = true;   % no MATLAB v. 5 legacy for ODE15I
end

if FcnHandlesUsed  % function handles used
  if isempty(tspan) || isempty(y0)
    error(message('MATLAB:odearguments:TspanOrY0NotSupplied', solver));
  end
  if length(tspan) < 2
    error(message('MATLAB:odearguments:SizeTspan', solver));
  end
  htspan = abs(tspan(2) - tspan(1));
  tspan = tspan(:);
  ntspan = length(tspan);
  t0 = tspan(1);
  next = 2;       % next entry in tspan
  tfinal = tspan(end);
  args = extras;                 % use f(t,y,p1,p2...)

else  % ode-file used   (ignored when solver == ODE15I)
  % Get default tspan and y0 from the function if none are specified.
  if isempty(tspan) || isempty(y0)
    if exist(ode)==2 && ( nargout(ode)<3 && nargout(ode)~=-1 )
      error(message('MATLAB:odearguments:NoDefaultParams', funstring( ode ), solver, funstring( ode )));
    end
    [def_tspan,def_y0,def_options] = feval(ode,[],[],'init',extras{:});
    if isempty(tspan)
      tspan = def_tspan;
    end
    if isempty(y0)
      y0 = def_y0;
    end
    options = odeset(def_options,options);
  end
  tspan = tspan(:);
  ntspan = length(tspan);
  if ntspan == 1    % Integrate from 0 to tspan
    t0 = 0;
    next = 1;       % Next entry in tspan.
  else
    t0 = tspan(1);
    next = 2;       % next entry in tspan
  end
  htspan = abs(tspan(next) - t0);
  tfinal = tspan(end);
  
  % The input arguments of f determine the args to use to evaluate f.
  if (exist(ode)==2)
    if (nargin(ode) == 2)
      args = {};                   % f(t,y)
    else
      args = [{''} extras];        % f(t,y,'',p1,p2...)
    end
  else  % MEX-files, etc.
    try
      args = [{''} extras];        % try f(t,y,'',p1,p2...)
      feval(ode,tspan(1),y0(:),args{:});
    catch
      args = {};                   % use f(t,y) only
    end
  end
end

y0 = y0(:);
neq = length(y0);

% Test that tspan is internally consistent.
if any(isnan(tspan))
  error(message('MATLAB:odearguments:TspanNaNValues'));
end
if t0 == tfinal
  error(message('MATLAB:odearguments:TspanEndpointsNotDistinct'));
end
tdir = sign(tfinal - t0);
if any( tdir*diff(tspan) <= 0 )
  error(message('MATLAB:odearguments:TspanNotMonotonic'));
end

f0 = feval(ode,t0,y0,args{:});   % ODE15I sets args{1} to yp0.
[m,n] = size(f0);
if n > 1
  error(message('MATLAB:odearguments:FoMustReturnCol', funstring( ode )));
elseif m ~= neq
    error(message('MATLAB:odearguments:SizeIC', funstring( ode ), m, neq, funstring( ode )));
end

% Determine the dominant data type
classT0 = class(t0);
classY0 = class(y0);
classF0 = class(f0);
if strcmp(solver,'ode15i')
  classYP0 = class(args{1});  % ODE15I sets args{1} to yp0.
  dataType = superiorfloat(t0,y0,args{1},f0);

  if ~( strcmp(classT0,dataType) && strcmp(classY0,dataType) && ...
        strcmp(classF0,dataType) && strcmp(classYP0,dataType))
    input1 = '''t0'', ''y0'', ''yp0''';
    input2 = '''f(t0,y0,yp0)''';
    warning(message('MATLAB:odearguments:InconsistentDataType',input1,input2,solver));
  end
else
  dataType = superiorfloat(t0,y0,f0);
  
  if ~( strcmp(classT0,dataType) && strcmp(classY0,dataType) && ...
        strcmp(classF0,dataType))
    input1 = '''t0'', ''y0''';
    input2 = '''f(t0,y0)''';
    warning(message('MATLAB:odearguments:InconsistentDataType',input1,input2,solver));
  end
end

% Get the error control options, and set defaults.
rtol = odeget(options,'RelTol',1e-3,'fast');
if (length(rtol) ~= 1) || (rtol <= 0)
  error(message('MATLAB:odearguments:RelTolNotPosScalar'));
end
if rtol < 100 * eps(dataType)
  rtol = 100 * eps(dataType);
  warning(message('MATLAB:odearguments:RelTolIncrease', sprintf( '%g', rtol )))
end
atol = odeget(options,'AbsTol',1e-6,'fast');
if any(atol <= 0)
  error(message('MATLAB:odearguments:AbsTolNotPos'));
end
normcontrol = strcmp(odeget(options,'NormControl','off','fast'),'on');
if normcontrol
  if length(atol) ~= 1
    error(message('MATLAB:odearguments:NonScalarAbsTol'));
  end
  normy = norm(y0);
else
  if (length(atol) ~= 1) && (length(atol) ~= neq)
    error(message('MATLAB:odearguments:SizeAbsTol', funstring( ode ), neq));
  end
  atol = atol(:);
  normy = [];
end
threshold = atol / rtol;

% By default, hmax is 1/10 of the interval.
safehmax = 16.0*eps(dataType)*max(abs(t0),abs(tfinal));  % 'inf' for tfinal = inf
defaulthmax = max(0.1*(abs(tfinal-t0)), safehmax);
hmax = min(abs(tfinal-t0), abs(odeget(options,'MaxStep',defaulthmax,'fast')));
if hmax <= 0
  error(message('MATLAB:odearguments:MaxStepLEzero'));
end
htry = abs(odeget(options,'InitialStep',[],'fast'));
if ~isempty(htry) && (htry <= 0)
  error(message('MATLAB:odearguments:InitialStepLEzero'));
end

odeFcn = ode;
end

function [haveeventfun,eventFcn,eventArgs,eventValue,teout,yeout,ieout] =...
    odeevents(FcnHandlesUsed,ode,t0,y0,options,extras)
%ODEEVENTS  Helper function for the events function in ODE solvers
%    ODEEVENTS initializes eventFcn to the events function, and creates a
%    cell-array of its extra input arguments. ODEEVENTS evaluates the events
%    function at(t0,y0).
%
%   See also ODE113, ODE15S, ODE23, ODE23S, ODE23T, ODE23TB, ODE45.

%   Jacek Kierzenka
%   Copyright 1984-2010 The MathWorks, Inc.

haveeventfun = 0;   % no Events function
eventArgs = [];
eventValue = [];
teout = [];
yeout = [];
ieout = [];

eventFcn = odeget(options,'Events',[],'fast');
if isempty(eventFcn)
  return
end

if FcnHandlesUsed     % function handles used
  haveeventfun = 1;   % there is an Events function
  eventArgs = extras;
  eventValue = feval(eventFcn,t0,y0,eventArgs{:});

else   % ode-file used
  switch lower(eventFcn)
    case 'on'
      haveeventfun = 1;   % there is an Events function
      eventFcn = ode;            % call ode(t,y,'events',p1,p2...)
      eventArgs = [{'events'}, extras];
      eventValue = feval(eventFcn,t0,y0,eventArgs{:});
    case 'off'
    otherwise
      error(message('MATLAB:odeevents:MustSetOnOrOff'))
  end
  
end
end

function newFcn = odefcncleanup(FcnHandleUsed,oldFcn,inputArgs)
% This helper function incorporates any input arguments for an ode function
% and creates a function handle from char inputs

%   Copyright 2017 The MathWorks, Inc.

if FcnHandleUsed
    if isempty(inputArgs)
        newFcn = oldFcn;
    else
        newFcn = @(t,y) oldFcn(t,y,inputArgs{:});
    end
else
    % Turn oldFcn string into function handle to avoid fevals
    [~,oldFcnFun] = evalc(['@' oldFcn]);
    if isempty(inputArgs)
        newFcn = @(t,y) oldFcnFun(t,y);
    else
        newFcn = @(t,y) oldFcnFun(t,y,inputArgs{:});
    end
end
end

function solver_output = odefinalize(solver, sol,...
                                     outfun, outargs,...
                                     printstats, statvect,...
                                     nout, tout, yout,...
                                     haveeventfun, teout, yeout, ieout,...
                                     interp_data)
%ODEFINALIZE Helper function called by ODE solvers at the end of integration.
%
%   See also ODE113, ODE15I, ODE15S, ODE23, ODE23S,
%            ODE23T, ODE23TB, ODE45, DDE23, DDESD.

%   Jacek Kierzenka
%   Copyright 1984-2005 The MathWorks, Inc.

if ~isempty(outfun)
  feval(outfun,[],[],'done',outargs{:});
end

% Return more stats for implicit solvers: ODE15i, ODE15s, ODE23s, ODE23t, ODE23tb
fullstats = (length(statvect) > 3);  % faster than 'switch' or 'ismember'

stats = struct('nsteps',statvect(1),'nfailed',statvect(2),'nfevals',statvect(3));
if fullstats
  stats.npds     = statvect(4);
  stats.ndecomps = statvect(5);
  stats.nsolves  = statvect(6);
else
  statvect(4:6) = 0;   % Backwards compatibility
end

if printstats
  fprintf(getString(message('MATLAB:odefinalize:LogSuccessfulSteps', sprintf('%g',stats.nsteps))));
  fprintf(getString(message('MATLAB:odefinalize:LogFailedAttempts', sprintf('%g',stats.nfailed))));
  fprintf(getString(message('MATLAB:odefinalize:LogFunctionEvaluations', sprintf('%g',stats.nfevals))));
  if fullstats
    fprintf(getString(message('MATLAB:odefinalize:LogPartialDerivatives', sprintf('%g',stats.npds))));
    fprintf(getString(message('MATLAB:odefinalize:LogLUDecompositions', sprintf('%g',stats.ndecomps))));
    fprintf(getString(message('MATLAB:odefinalize:LogSolutionsOfLinearSystems', sprintf('%g',stats.nsolves))));
  end
end

solver_output = {};

if (nout > 0) % produce output
  if isempty(sol) % output [t,y,...]
    solver_output{1} = tout(1:nout).';
    solver_output{2} = yout(:,1:nout).';
    if haveeventfun
      solver_output{3} = teout.';
      solver_output{4} = yeout.';
      solver_output{5} = ieout.';
    end
    solver_output{end+1} = statvect(:);  % Column vector
  else % output sol
    % Add remaining fields
    sol.x = tout(1:nout);
    sol.y = yout(:,1:nout);
    if haveeventfun
      sol.xe = teout;
      sol.ye = yeout;
      sol.ie = ieout;
    end
    sol.stats = stats;
    switch solver
     case {'dde23','ddesd'}
      [history,ypout] = deal(interp_data{:});
      sol.yp = ypout(:,1:nout);
      if isstruct(history)
        sol.x = [history.x sol.x];
        sol.y = [history.y sol.y];
        sol.yp = [history.yp sol.yp];
        if isfield(history,'xe')
          if isfield(sol,'xe')
            sol.xe = [history.xe sol.xe];
            sol.ye = [history.ye sol.ye];
            sol.ie = [history.ie sol.ie];
          else
            sol.xe = history.xe;
            sol.ye = history.ye;
            sol.ie = history.ie;
          end
        end
      end
     case 'ode45'
      [f3d,idxNonNegative] = deal(interp_data{:});
      sol.idata.f3d = f3d(:,:,1:nout);
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode15s'
      [kvec,dif3d,idxNonNegative] = deal(interp_data{:});
      sol.idata.kvec = kvec(1:nout);
      maxkvec = max(sol.idata.kvec);
      sol.idata.dif3d = dif3d(:,1:maxkvec+2,1:nout);
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode113'
      [klastvec,phi3d,psi2d,idxNonNegative] = deal(interp_data{:});
      sol.idata.klastvec = klastvec(1:nout);
      kmax = max(sol.idata.klastvec);
      sol.idata.phi3d = phi3d(:,1:kmax+1,1:nout);
      sol.idata.psi2d = psi2d(1:kmax,1:nout);
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode23'
      [f3d,idxNonNegative] = deal(interp_data{:});
      sol.idata.f3d = f3d(:,:,1:nout);
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode23s'
      [k1data,k2data] = deal(interp_data{:});
      sol.idata.k1 = k1data(:,1:nout);
      sol.idata.k2 = k2data(:,1:nout);
     case 'ode23t'
      [zdata,znewdata,idxNonNegative] = deal(interp_data{:});
      sol.idata.z = zdata(:,1:nout);
      sol.idata.znew = znewdata(:,1:nout);
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode23tb'
      [t2data,y2data,idxNonNegative] = deal(interp_data{:});
      sol.idata.t2 = t2data(1:nout);
      sol.idata.y2 = y2data(:,1:nout);
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode15i'
      [kvec,ypfinal] = deal(interp_data{:});
      sol.idata.kvec = kvec(1:nout);
      sol.extdata.ypfinal = ypfinal;
     otherwise
      error(message('MATLAB:odefinalize:UnrecognizedSolver', solver));
    end
    if strcmp(solver,'dde23') || strcmp(solver,'ddesd')
      solver_output = sol;
    else
      solver_output{1} = sol;
    end
  end
end
end

function [Jconstant,Jfcn,Jargs,Joptions] = ...
    odejacobian(fcnHandlesUsed,ode,t0,y0,options,extras)
%ODEJACOBIAN  Helper function for the Jacobian function in ODE solvers
%    ODEJACOBIAN determines whether the Jacobian is constant and if so,
%    returns its value as Jfcn. If an analytical Jacobian is available from
%    a function, ODEJACOBIAN initializes Jfcn and creates a cell array of
%    additional input arguments. For numerical Jacobian, ODEJACOBIAN tries to
%    extract JPattern and sets JOPTIONS for use with ODENUMJAC.
%
%   See also ODE15S, ODE23S, ODE23T, ODE23TB, ODENUMJAC.

%   Jacek Kierzenka
%   Copyright 1984-2009 The MathWorks, Inc.

Jconstant = strcmp(odeget(options,'JConstant','off','fast'),'on');
Jfcn = [];
Jargs = {};
Joptions = [];

Janalytic = false;

if fcnHandlesUsed
  Jfcn = odeget(options,'Jacobian',[],'fast');
  if ~isempty(Jfcn)
    if isnumeric(Jfcn)
      Jconstant = true;
    else
      Janalytic = true;
      Jargs = extras;
    end
  end
else  % ode-file used
  joption = odeget(options,'Jacobian','off','fast');
  switch lower(joption)
    case 'on'    % ode(t,y,'jacobian',p1,p2...)
      Janalytic = true;
      Jfcn = ode;
      Jargs = [{'jacobian'} extras];
    case 'off'   % use odenumjac
    otherwise
      error(message('MATLAB:odejacobian:InvalidJOption', joption));
  end
end

if ~Janalytic   % odenumjac will be used
  Joptions.diffvar  = 2;       % df(t,y)/dy
  Joptions.vectvars = [];
  vectorized = strcmp(odeget(options,'Vectorized','off','fast'),'on');
  if vectorized
    Joptions.vectvars = 2;     % f(t,[y1,y2]) = [f(t,y1), f(t,y2)]
  end
  
  atol = odeget(options,'AbsTol',1e-6,'fast');
  Joptions.thresh = zeros(size(y0))+ atol(:);
  Joptions.fac  = [];
  
  if fcnHandlesUsed
    jpattern = odeget(options,'JPattern',[],'fast');
  else  % ode-file used
    jp_option = odeget(options,'JPattern','off','fast');
    switch lower(jp_option)
      case 'on'
        jpattern = feval(ode,[],[],'jpattern',extras{:});
      case 'off'  % no pattern provided
        jpattern = [];
      otherwise
        error(message('MATLAB:odejacobian:InvalidJpOption', jp_option));
    end
  end
  if ~isempty(jpattern)
    Joptions.pattern = jpattern;
    Joptions.g = colgroup(jpattern);
  end
end
end
    
function [massType, massM, massFcn, massArgs, dMoptions] = ...
    odemass(FcnHandlesUsed,ode,t0,y0,options,extras)
%ODEMASS  Helper function for the mass matrix function in ODE solvers
%    ODEMASS determines the type of the mass matrix, initializes massFcn to
%    the mass matrix function and creates a cell-array of extra input
%    arguments. ODEMASS evaluates the mass matrix at(t0,y0).
%
%   See also ODE113, ODE15S, ODE23, ODE23S, ODE23T, ODE23TB, ODE45.

%   Jacek Kierzenka
%   Copyright 1984-2011 The MathWorks, Inc.

massType = 0;
massFcn = [];
massArgs = {};
massM = speye(length(y0));
dMoptions = [];    % options for odenumjac computing d(M(t,y)*v)/dy

if FcnHandlesUsed     % function handles used
  Moption = odeget(options,'Mass',[],'fast');
  if isempty(Moption)
    return    % massType = 0
  elseif isnumeric(Moption)
    massType = 1;
    massM = Moption;
  else % try feval
    massFcn = Moption;
    massArgs = extras;
    Mstdep = odeget(options,'MStateDependence','weak','fast');
    switch lower(Mstdep)
      case 'none'
        massType = 2;
      case 'weak'
        massType = 3;
      case 'strong'
        massType = 4;
        
        dMoptions.diffvar  = 3;       % d(odeMxV(Mfun,t,y)/dy
        dMoptions.vectvars = [];
        
        atol = odeget(options,'AbsTol',1e-6,'fast');
        dMoptions.thresh = zeros(size(y0))+ atol(:);
        
        dMoptions.fac  = [];
        
        Mvs = odeget(options,'MvPattern',[],'fast');
        if ~isempty(Mvs)
          dMoptions.pattern = Mvs;
          dMoptions.g = colgroup(Mvs);
        end
                  
      otherwise
        error(message('MATLAB:odemass:MStateDependenceMassType'));
    end
    if massType > 2   % state-dependent
      massM = feval(massFcn,t0,y0,massArgs{:});
    else   % time-dependent only
      massM = feval(massFcn,t0,massArgs{:});
    end
  end
  
else % ode-file
  mass = lower(odeget(options,'Mass','none','fast'));

  switch(mass)
    case 'none', return;  % massType = 0
    case 'm', massType = 1;
    case 'm(t)', massType = 2;
    case 'm(t,y)', massType = 3;
    otherwise
      error(message('MATLAB:odemass:InvalidMassProp', mass));
  end
  massFcn = ode;
  massArgs = [{'mass'}, extras];
  massM = feval(massFcn,t0,y0,massArgs{:});

end
end

function [odeFcn,odeArgs] = odemassexplicit( FcnHandlesUsed,massType,odeFcn,...
                                             odeArgs,massFcn,massM)
%ODEMASSEXPLICIT  Helper function for handling the mass matrix
%   For explicit ODE solvers -- incorporate the mass matrix into the ODE
%   function.
%
%   See also ODE113, ODE23, ODE45.

%   Jacek Kierzenka
%   Copyright 1984-2006 The MathWorks, Inc.

if FcnHandlesUsed
    switch massType
      case 1  % use LU factors of constant M
        if issparse(massM)
            [massL,massU,massP,massQ,massR] = lu(massM);
            odeArgs = [{odeFcn,massL,massU,massP,massQ,massR},odeArgs];
            odeFcn = @ExplicitSolverHandleMass1sparse;
        else % M full
            [massL,massU,massp] = lu(massM,'vector');
            odeArgs = [{odeFcn,massL,massU,massp},odeArgs];
            odeFcn = @ExplicitSolverHandleMass1;
        end
      case 2
        odeArgs = [{odeFcn,massFcn},odeArgs];
        odeFcn = @ExplicitSolverHandleMass2;
      otherwise % case {3,4}
        odeArgs = [{odeFcn,massFcn},odeArgs];
        odeFcn = @ExplicitSolverHandleMass34;
    end
else % ode-file:  F(t,y,'mass',p1,p2...)
    if massType == 1   % use LU factors of constant M
        if issparse(massM)
            [massL,massU,massP,massQ,massR] = lu(massM);
            odeArgs = [{odeFcn,massL,massU,massP,massQ,massR},odeArgs];
            odeFcn = @ExplicitSolverHandleMass1sparse;
        else % M full
            [massL,massU,massp] = lu(massM,'vector');
            odeArgs = [{odeFcn,massL,massU,massp},odeArgs];
            odeFcn = @ExplicitSolverHandleMass1;
        end
    else
        odeArgs = [{odeFcn},odeArgs];
        odeFcn = @ExplicitSolverHandleMassOld;
    end
end
end

% --------------------------------------------------------------------------

function yp = ExplicitSolverHandleMass1(t,y,odeFcn,L,U,p,varargin)
  ode = feval(odeFcn,t,y,varargin{:});
  yp = U \ (L \ ode(p));
end

% --------------------------------------------------------------------------

function yp = ExplicitSolverHandleMass1sparse(t,y,odeFcn,L,U,P,Q,R,varargin)
  yp = Q *( U \ (L \ (P * (R \ feval(odeFcn,t,y,varargin{:})))));
end
 
% --------------------------------------------------------------------------
  
function yp = ExplicitSolverHandleMass2(t,y,odeFcn,massFcn,varargin)
  yp = feval(massFcn,t,varargin{:}) \ feval(odeFcn,t,y,varargin{:});
end
  
% --------------------------------------------------------------------------

function yp = ExplicitSolverHandleMass34(t,y,odeFcn,massFcn,varargin)
  yp = feval(massFcn,t,y,varargin{:}) \ feval(odeFcn,t,y,varargin{:});
end

% --------------------------------------------------------------------------
  
function yp = ExplicitSolverHandleMassOld(t,y,odeFcn,varargin)
  yp = feval(odeFcn,t,y,'mass',varargin{2:end}) \ ...
       feval(odeFcn,t,y,varargin{:});
end
  
% --------------------------------------------------------------------------

function out = odemxv(Mfun,t,y,v,varargin)
%ODEMXV  Helper function -- evaluates Mfun(t,y)*v
%   Used to get d(M(t,y)*v)/dy when the property MStateDependence is 'strong'
%
%   See also DAEIC3, ODE15S, ODE23T, ODE23TB.

%   Jacek Kierzenka, Lawrence Shampine
%   Copyright 1984-2002 The MathWorks, Inc.

out = feval(Mfun,t,y,varargin{:})*v;
end

function [odeFcn,thresholdNonNegative] = odenonnegative(ode,y0,threshold,idxNonNegative)
%ODENONNEGATIVE  Helper function for handling nonnegative solution constraints
%   Modify the derivative function to prevent the solution from crossing zero.
%
%   See also ODE113, ODE15S, ODE23, ODE23T, ODE23TB, ODE45.

%   Jacek Kierzenka
%   Copyright 1984-2010 The MathWorks, Inc.

neq = numel(y0);
thresholdNonNegative = [];
if any( (idxNonNegative < 1) | (idxNonNegative > neq) )
  error(message('MATLAB:odenonnegative:NonNegativeIndicesInvalid'));
end
if any(y0(idxNonNegative) < 0)
  error(message('MATLAB:odenonnegative:NonNegativeViolatedAtT0'));
end
if length(threshold) == 1
  thresholdNonNegative = threshold(ones(size(idxNonNegative)));
else
  thresholdNonNegative = threshold(idxNonNegative);
end
thresholdNonNegative = thresholdNonNegative(:);
odeFcn = @local_odeFcn_nonnegative;

% -----------------------------------------------------------
% Nested function: ODE with nonnegativity constraints imposed
%
  function yp = local_odeFcn_nonnegative(t,y,varargin)
    yp = feval(ode,t,y,varargin{:});
    ndx = idxNonNegative( find(y(idxNonNegative) <= 0) );
    yp(ndx) = max(yp(ndx),0);
  end  % local_odeFcn_nonnegative
% -----------------------------------------------------------

end  % odenonnegative

function [dFdy,fac,nfevals,nfcalls] = odenumjac(F,Fargs,Fvalue,options)
%ODENUMJAC Numerically compute the Jacobian dF/dY of function F(...,Y,...).
%   [DFDY,FAC] = ODENUMJAC(F,FARGS,FVALUE,OPTIONS) numerically computes
%   the Jacobian of function F with respect to variable Y, returning it
%   as a matrix DFDY. F could be a function of several variables. It must
%   return a column vector. The arguments of F are specified in a cell
%   array FARGS. Vector FVALUE contains F(FARGS{:}).
%   The structure OPTIONS must have the following fields: DIFFVAR, VECTVARS,
%   THRESH, and FAC. For sparse Jacobians, OPTIONS must also have fields
%   PATTERNS and G. The field OPTIONS.DIFFVAR is the index of the
%   differentiation variable, Y = FARGS{DIFFVAR}. For a function F(t,x),
%   set DIFFVAR to 1 to compute DF/Dt, or to 2 to compute DF/Dx.
%   ODENUMJAC takes advantage of vectorization, i.e., when several values F
%   can be obtained with one function evaluation. Set OPTIONS.VECTVAR
%   to the indices of vectorized arguments: VECTVAR = [2] indicates that
%   F(t,[x1 y2 ...]) returns [F(t,x1) F(t,x2) ...], while VECTVAR = [1,2]
%   indicates that F([t1 t2 ...],[x1 x2 ...]) returns [F(t1,x1) F(t2,x2) ...].
%   OPTIONS.THRESH provides a threshold of significance for Y, i.e.
%   the exact value of a component Y(i) with abs(Y(i)) < THRESH(i) is not
%   important. All components of THRESH must be positive. Column FAC is
%   working storage. On the first call, set OPTIONS.FAC to []. Do not alter
%   the returned value between calls.
%
%   [DFDY,FAC] = ODENUMJAC(F,FARGS,FVALUE,OPTIONS) computes a sparse matrix
%   DFDY if the fields OPTIONS.PATTERN and OPTIONS.G are present.
%   PATTERN is a non-empty sparse matrix of 0's and 1's. A value of 0 for
%   PATTERN(i,j) means that component i of the function F(...,Y,...) does not
%   depend on component j of vector Y (hence DFDY(i,j) = 0).  Column vector
%   OPTIONS.G is an efficient column grouping, as determined by COLGROUP(PATTERN).
%
%   [DFDY,FAC,NFEVALS,NFCALLS] = ODENUMJAC(...) returns the number of values
%   F(FARGS{:}) computed while forming dFdY (NFEVALS) and the number of calls
%   to the function F (NFCALLS). If F is not vectorized, NFCALLS equals NFEVALS.
%
%   Although ODENUMJAC was developed specifically for the approximation of
%   partial derivatives when integrating a system of ODE's, it can be used
%   for other applications.  In particular, when the length of the vector
%   returned by F(...,Y,...) is different from the length of Y, DFDY is
%   rectangular.
%
%   See also COLGROUP.

%   ODENUMJAC is an implementation of an exceptionally robust scheme due to
%   Salane for the approximation of partial derivatives when integrating
%   a system of ODEs, Y' = F(T,Y). It is called when the ODE code has an
%   approximation Y at time T and is about to step to T+H.  The ODE code
%   controls the error in Y to be less than the absolute error tolerance
%   ATOL = THRESH.  Experience computing partial derivatives at previous
%   steps is recorded in FAC.  A sparse Jacobian is computed efficiently
%   by using COLGROUP(S) to find groups of columns of DFDY that can be
%   approximated with a single call to function F.  COLGROUP tries two
%   schemes (first-fit and first-fit after reverse COLAMD ordering) and
%   returns the better grouping.
%
%   D.E. Salane, "Adaptive Routines for Forming Jacobians Numerically",
%   SAND86-1319, Sandia National Laboratories, 1986.
%
%   T.F. Coleman, B.S. Garbow, and J.J. More, Software for estimating
%   sparse Jacobian matrices, ACM Trans. Math. Software, 11(1984)
%   329-345.
%
%   L.F. Shampine and M.W. Reichelt, The MATLAB ODE Suite, SIAM Journal on
%   Scientific Computing, 18-1, 1997.

%   Mark W. Reichelt and Lawrence F. Shampine, 3-28-94
%   Copyright 1984-2012 The MathWorks, Inc.

% Options
diffvar = options.diffvar;
vectvar = options.vectvars;
thresh  = options.thresh;
fac     = options.fac;

% Full or sparse Jacobian.
fullJacobian = true;
if isfield(options,'pattern')
  fullJacobian = false;
  S = options.pattern;
  g = options.g;
end
  
% The differentiation variable.
y  = Fargs{diffvar};

% Initialize.
Fvalue = Fvalue(:);
classF = class(Fvalue);
br = eps(classF) ^ (0.875);
bl = eps(classF) ^ (0.75);
bu = eps(classF) ^ (0.25);
classY = class(y);
facmin = eps(classY) ^ (0.78);
facmax = 0.1;
ny = length(y);
nF = length(Fvalue);
if isempty(fac)
  fac = sqrt(eps(classY)) + zeros(ny,1,classY);
end

% Select an increment del for a difference approximation to
% column j of dFdy.  The vector fac accounts for experience
% gained in previous calls to numjac.
yscale = max(abs(y),thresh);
del = (y + fac .* yscale) - y;
for j = find(del == 0)'
  while true
    if fac(j) < facmax
      fac(j) = min(100*fac(j),facmax);
      del(j) = (y(j) + fac(j)*yscale(j)) - y(j);
      if del(j)
        break
      end
    else
      del(j) = thresh(j);
      break;
    end
  end
end
if nF == ny
  s = (sign(Fvalue) >= 0);
  del = (s - (~s)) .* abs(del);         % keep del pointing into region
end

% Form a difference approximation to all columns of dFdy.
if fullJacobian                           % generate full matrix dFdy
  ydel = y(:,ones(1,ny)) + diag(del);
  if isempty(vectvar)
    % non-vectorized
    Fdel = zeros(nF,ny);
    for j = 1:ny
      Fdel(:,j) = feval(F,Fargs{1:diffvar-1},ydel(:,j),Fargs{diffvar+1:end});
    end
    nfcalls = ny;                       % stats
  else
    % Expand arguments.  Need to preserve the original (non-expanded)
    % Fargs in case of correcting columns (done one column at a time).
    Fargs_expanded = Fargs;
    Fargs_expanded{diffvar} = ydel;
    vectvar = setdiff(vectvar,diffvar);
    for i=1:length(vectvar)
      Fargs_expanded{vectvar(i)} = repmat(Fargs{vectvar(i)},1,ny);
    end
    Fdel = feval(F,Fargs_expanded{:});
    nfcalls = 1;                        % stats
  end
  nfevals = ny;                         % stats (at least one per loop)
  Fdiff = Fdel - Fvalue(:,ones(1,ny));
  dFdy = Fdiff * diag(1 ./ del);
  [Difmax,Rowmax] = max(abs(Fdiff),[],1);
  % If Fdel is a column vector, then index is a scalar, so indexing is okay.
  absFdelRm = abs(Fdel((0:ny-1)*nF + Rowmax));

else                    % sparse dFdy with structure S and column grouping g
  nzcols = find(g > 0); % g==0 for all-zero columns in sparsity pattern

  if isempty(nzcols)  % No columns requested -- early exit
    dFdy = sparse([],[],[],nF,ny);
    nfcalls = 0;                        % stats
    nfevals = 0;                        % stats
    return;
  end

  ng = max(g);
  one2ny = (1:ny)';
  ydel = y(:,ones(1,ng));
  
  i = (g(nzcols)-1)*ny + nzcols;
  ydel(i) = ydel(i) + del(nzcols);
  if isempty(vectvar)
    % non-vectorized
    Fdel = zeros(nF,ng);
    for j = 1:ng
      Fdel(:,j) = feval(F,Fargs{1:diffvar-1},ydel(:,j),Fargs{diffvar+1:end});
    end
    nfcalls = ng;                       % stats
  else
    % Expand arguments.  Need to preserve the original (non-expanded)
    % Fargs in case of correcting columns (done one column at a time).
    Fargs_expanded = Fargs;
    Fargs_expanded{diffvar} = ydel;
    vectvar = setdiff(vectvar,diffvar);
    for i=1:length(vectvar)
      Fargs_expanded{vectvar(i)} = repmat(Fargs{vectvar(i)},1,ng);
    end
    Fdel = feval(F,Fargs_expanded{:});
    nfcalls = 1;                        % stats
  end
  nfevals = ng;                         % stats (at least one per column)
  Fdiff = Fdel - Fvalue(:,ones(1,ng));
  [i, j] = find(S);
  if ~isempty(i)
    i = i(:); % ensure that i is a column vector (S could be a row vector)
  end
  Fdiff = sparse(i,j,Fdiff((g(j)-1)*nF + i),nF,ny);
  dFdy = Fdiff * sparse(one2ny,one2ny,1 ./ del,ny,ny);
  [Difmax,Rowmax] = max(abs(Fdiff),[],1);
  Difmax = full(Difmax);
  % If ng==1, then Fdel is a column vector although index may be a row vector.
  absFdelRm = zeros(1,ny);
  absFdelRm(nzcols) = abs(Fdel((g(nzcols)-1)*nF + Rowmax(nzcols)'));
end

% Adjust fac for next call to numjac.
absFvalue = abs(Fvalue);
absFvalueRm = absFvalue(Rowmax);              % not a col vec if absFvalue scalar
absFvalueRm = absFvalueRm(:)';                % ensure that absFvalueRm is a row vector
absFdelRm = absFdelRm(:)';                    % ensure that absFdelRm is a row vector
j = ((absFdelRm ~= 0) & (absFvalueRm ~= 0)) | (Difmax == 0);

if ~fullJacobian
    j = j & (g(:)'>0);   % refine only requested columns
end
    
if any(j)
  ydel = y;
  Fscale = max(absFdelRm,absFvalueRm);

  % If the difference in f values is so small that the column might be just
  % roundoff error, try a bigger increment.
  k1 = (Difmax <= br*Fscale);           % Difmax and Fscale might be zero
  for k = find(j & k1)
    tmpfac = min(sqrt(fac(k)),facmax);
    del = (y(k) + tmpfac*yscale(k)) - y(k);
    if (tmpfac ~= fac(k)) && (del ~= 0)
      if nF == ny
        if Fvalue(k) >= 0                  % keep del pointing into region
          del = abs(del);
        else
          del = -abs(del);
        end
      end
        
      ydel(k) = y(k) + del;
      Fargs{diffvar} = ydel;
      fdel = feval(F,Fargs{:});
      nfevals = nfevals + 1;            % stats
      nfcalls = nfcalls + 1;            % stats
      ydel(k) = y(k);
      fdiff = fdel(:) - Fvalue;
      tmp = fdiff ./ del;
      
      [difmax,rowmax] = max(abs(fdiff));
      if tmpfac * norm(tmp,inf) >= norm(dFdy(:,k),inf);
        % The new difference is more significant, so
        % use the column computed with this increment.
        if fullJacobian
          dFdy(:,k) = tmp;
        else
          i = find(S(:,k));
          if ~isempty(i)
            dFdy(i,k) = tmp(i);
          end
        end
  
        % Adjust fac for the next call to numjac.
        fscale = max(abs(fdel(rowmax)),absFvalue(rowmax));
          
        if difmax <= bl*fscale
          % The difference is small, so increase the increment.
          fac(k) = min(10*tmpfac, facmax);
        elseif difmax > bu*fscale
          % The difference is large, so reduce the increment.
          fac(k) = max(0.1*tmpfac, facmin);
        else
          fac(k) = tmpfac;
        end
      end
    end
  end
  
  % If the difference is small, increase the increment.
  k = find(j & ~k1 & (Difmax <= bl*Fscale));
  if ~isempty(k)
    fac(k) = min(10*fac(k), facmax);
  end

  % If the difference is large, reduce the increment.
  k = find(j & (Difmax > bu*Fscale));
  if ~isempty(k)
    fac(k) = max(0.1*fac(k), facmin);
  end
end
end

function [tout,yout,iout,vnew,stop] = ...
    odezero(ntrpfun,eventfun,eventargs,v,t,y,tnew,ynew,t0,varargin)
%ODEZERO Locate any zero-crossings of event functions in a time step.
%   ODEZERO is an event location helper function for the ODE Suite.  ODEZERO
%   uses Regula Falsi and information passed from the ODE solver to locate
%   any zeros in the half open time interval (T,TNEW] of the event functions
%   coded in eventfun.
%
%   See also ODE45, ODE23, ODE113, ODE15S, ODE23S, ODE23T, ODE23TB.

%   Mark W. Reichelt, Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2019 The MathWorks, Inc.

% Initialize.
tol = 128*max(eps(t),eps(tnew));
tol = min(tol, abs(tnew - t));
tout = [];
yout = [];
iout = [];
tdir = sign(tnew - t);
stop = 0;
rmin = realmin;

% Set up tL, tR, yL, yR, vL, vR, isterminal and direction.
tL = t;
yL = y;
vL = v;
[vnew,isterminal,direction] = feval(eventfun,tnew,ynew,eventargs{:});
if isempty(direction)
  direction = zeros(size(vnew));   % zeros crossings in any direction
end
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
    indzc = find((sign(vL) ~= sign(vR)) & (direction .* (vR - vL) >= 0));
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
    ytry = feval(ntrpfun,ttry,t,y,tnew,ynew,varargin{:});
    vtry = feval(eventfun,ttry,ytry,eventargs{:});

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
  end

  j = ones(1,length(indzc));
  tout = [tout, tR(j)];
  yout = [yout, yR(:,j)];
  iout = [iout, indzc(:)'];
  if any(isterminal(indzc))
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
    yL = feval(ntrpfun,tL,t,y,tnew,ynew,varargin{:});
    vL = feval(eventfun,tL,yL,eventargs{:});
    tR = tnew; yR = ynew; vR = vnew;
  end
end
end

function [yinterp,ypinterp] = ntrp45(tinterp,t,y,~,~,h,f,idxNonNegative)
%NTRP45  Interpolation helper function for ODE45.
%   YINTERP = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F,IDX) uses data computed in ODE45
%   to approximate the solution at time TINTERP.  TINTERP may be a scalar
%   or a row vector.
%   The arguments TNEW and YNEW do not affect the computations. They are
%   required for consistency of syntax with other interpolation functions.
%   Any values entered for TNEW and YNEW are ignored.
%
%   [YINTERP,YPINTERP] = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F,IDX) returns also the
%   derivative of the polynomial approximating the solution.
%
%   IDX has indices of solution components that must be non-negative. Negative
%   YINTERP(IDX) are replaced with zeros and the derivative YPINTERP(IDX) is
%   set to zero.
%
%   See also ODE45, DEVAL.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2009 The MathWorks, Inc.

BI = [
    1       -183/64      37/12       -145/128
    0          0           0            0
    0       1500/371    -1000/159    1000/371
    0       -125/32       125/12     -375/64
    0       9477/3392   -729/106    25515/6784
    0        -11/7        11/3        -55/28
    0         3/2         -4            5/2
    ];

s = (tinterp - t)/h;
yinterp = y(:,ones(size(tinterp))) + f*(h*BI)*cumprod([s;s;s;s]);

ypinterp = [];
if nargout > 1
  ypinterp = f*BI*[ ones(size(s)); cumprod([2*s;3/2*s;4/3*s])];
end

% Non-negative solution
if ~isempty(idxNonNegative)
  idx = find(yinterp(idxNonNegative,:)<0); % vectorized
  if ~isempty(idx)
    w = yinterp(idxNonNegative,:);
    w(idx) = 0;
    yinterp(idxNonNegative,:) = w;
    if nargout > 1   % the derivative
      w = ypinterp(idxNonNegative,:);
      w(idx) = 0;
      ypinterp(idxNonNegative,:) = w;
    end
  end
end
end

function yinterp = ntrp45split(tinterp,t,y,h,f1,f3,f4,f5,f6,f7,idxNonNegative)
%NTRP45SPLIT  Interpolation helper function for ODE45.
%   YINTERP = NTRP45SPLIT(TINTERP,T,Y,H,F1,F3,F4,F5,F6,F7,IDX) uses data
%   computed in ODE45 to approximate the solution at time TINTERP.  TINTERP
%   may be a scalar or a row vector.
%
%   IDX has indices of solution components that must be non-negative.
%   Negative YINTERP(IDX) are replaced with zeros.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2017 The MathWorks, Inc.

% Define constants as scalars
bi12 = -183/64;   bi13 = 37/12;     bi14 = -145/128;
bi32 = 1500/371;  bi33 = -1000/159; bi34 = 1000/371;
bi42 = -125/32;   bi43 = 125/12;    bi44 = -375/64;
bi52 = 9477/3392; bi53 = -729/106;  bi54 = 25515/6784;
bi62 = -11/7;     bi63 = 11/3;      bi64 = -55/28;
bi72 = 3/2;       bi73 = -4;        bi74 = 5/2;

s = (tinterp - t)/h;

% Preallocate array then use for loop to iterate
yinterp = zeros(size(y, 1), size(s, 2));
for jj=1:size(s, 2)
    sj = s(jj);
    sj2 = sj.*sj;
    bs1 = (sj + sj2.*(bi12 + sj.*(bi13 + bi14*sj)));
    bs3 = (     sj2.*(bi32 + sj.*(bi33 + bi34*sj)));
    bs4 = (     sj2.*(bi42 + sj.*(bi43 + bi44*sj)));
    bs5 = (     sj2.*(bi52 + sj.*(bi53 + bi54*sj)));
    bs6 = (     sj2.*(bi62 + sj.*(bi63 + bi64*sj)));
    bs7 = (     sj2.*(bi72 + sj.*(bi73 + bi74*sj)));
    
    yinterp(:,jj) = y + h*(f1.*bs1 + f3.*bs3 + f4.*bs4 + f5.*bs5 + f6.*bs6 + f7.*bs7);
end

% Non-negative solution
if ~isempty(idxNonNegative)
  idx = find(yinterp(idxNonNegative,:)<0); % vectorized
  if ~isempty(idx)
    w = yinterp(idxNonNegative,:);
    w(idx) = 0;
    yinterp(idxNonNegative,:) = w;
  end
end
end
