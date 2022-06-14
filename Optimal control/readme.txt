The MATLAB code files in this folder implement the methods for tractable global solutions to single-input optimal control problems in the presence of plant-model mismatch that are described in the papers "Efficient Global Solutions to Single-Input Optimal Control Problems via Approximation by Sum-of-Squares Polynomials", "Tractable Global Solutions to Chance-Constrained Bayesian Optimal Experiment Design for Arbitrary Prior and Noise Distributions", and "Data-Driven Adaptive Optimal Control Under Model Uncertainty: An Application to Cold Atmospheric Plasmas". Other code files are provided for further examples.

The use of the code in this folder assumes a correct installation of the SDP solver MOSEK (https://www.mosek.com) and the algorithmic differentiation toolbox CasADi (https://web.casadi.org), in addition to a running version of MATLAB with the MEX compiler, the Symbolic Toolbox, the Optimization Toolbox, and the Parallel Computing Toolbox.

Before the first use of the code files in this folder, start MATLAB in this folder and execute the following commands:
cd ..
cd 'multiv_pols_min'
make_mex
addpath('.')
cd ..
cd 'ode45_rp'
make_mex
addpath('.')
cd ..
cd 'Optimal control'
savepath

The code files provide the results for the aforementioned papers and further examples, as follows:
- The code file test_hessian_pyrrole_bc.m provides the results for the reaction system of pyrrole acetoacetylation in Section VI.A of the paper "Efficient Global Solutions to Single-Input Optimal Control Problems via Approximation by Sum-of-Squares Polynomials", namely Table I and Fig. 1. The code file test_hessian_goddard.m provides the results for the Goddard problem in Section VI.B of the paper "Efficient Global Solutions to Single-Input Optimal Control Problems via Approximation by Sum-of-Squares Polynomials", namely Table II and Fig. 2. The code files test_hessian_time.m and test_hessian_cost.m provide similar results for two additional examples with the chlorination system in a previous paper "Dynamic Optimization of Reaction Systems via Exact Parsimonious Input Parameterization". The code file test_casadi_pyrrole_bc.m provides results of an adaptive approach for the reaction system of pyrrole acetoacetylation in Section VI.A of the paper "Efficient Global Solutions to Single-Input Optimal Control Problems via Approximation by Sum-of-Squares Polynomials", as well as the results of model predictive control by using a parsimonious parameterization approach, single shooting, and multiple shooting.
- The code file test_hessian_lv_fim.m provides the results for the Lotka-Volterra system in Section 7.1 of the paper "Tractable Global Solutions to Chance-Constrained Bayesian Optimal Experiment Design for Arbitrary Prior and Noise Distributions", namely Table 1. The code file test_hessian_lv_mc.m provides the results for the Lotka-Volterra system in Section 7.2 of the paper "Tractable Global Solutions to Chance-Constrained Bayesian Optimal Experiment Design for Arbitrary Prior and Noise Distributions", namely Table 2. The code files test_hessian_lin_fim.m and test_hessian_lin_mc.m can be used to obtain similar results for a linear system.
- The code file test_hessian_appj.m provides the results that would be obtained if the true APPJ were perfectly known in Section IV of the paper "Data-Driven Adaptive Optimal Control Under Model Uncertainty: An Application to Cold Atmospheric Plasmas". The code file test_appj_track.m provides the results in the presence of plant-model mismatch and disturbances for the APPJ in Section IV of the paper "Data-Driven Adaptive Optimal Control Under Model Uncertainty: An Application to Cold Atmospheric Plasmas". The code file test_casadi_appj.m provides results of an adaptive approach for the APPJ in Section IV of the paper "Data-Driven Adaptive Optimal Control Under Model Uncertainty: An Application to Cold Atmospheric Plasmas", as well as the results of model predictive control by using a parsimonious parameterization approach, single shooting, and multiple shooting.
- The code file test_hessian_lv_lqr.m provides the solution to a problem related to the Lotka-Volterra system for an example provided with CasADi.
- The code file test_hessian_example.m provides the results for an example of an optimal control problem with a singular arc.
- The code file test_hessian_simple.m provides useful examples for testing the correct computation of first-order and second-order partial derivatives with respect to different variables in various scenarios.

To achieve global optimality in a tractable way, the proposed procedure (i) approximates the optimal control problem (OCP) as a set of polynomial optimization problems (POPs), (ii) computes the global solution to each POP by reformulating it as a semidefinite program (SDP) via the concept of sum-of-squares polynomials, and (iii) uses the solution to the POPs to compute a global solution to the OCP. The approximation of the OCP as a set of POPs is implemented by a library of MATLAB code files composed of the files hessian_seq.m, hessian_approx.m, hessian_eval.m, hessian_svm.m, and hessian_regress.m. The computation of the global solution to each POP by reformulating it as an SDP via the concept of sum-of-squares polynomials is implemented by a library of MATLAB code files composed of the file hessian_optim.m, as well as the files multiv_pols_min.m, sdp_mosek.m, calc_mon.m, solve_nonlin.m, and auxiliary C files in the folder multiv_pols_min. These libraries are regarded as black boxes, thus no implementation details are given here for the sake of simplicity. The use of the solution to the POPs to compute a global solution to the OCP is implemented by the file hessian_fmincon.m. The code files for (i) and (iii) use the file hessian_calc.m, as well as the file ode45_rp.m and auxiliary C files in the folder ode45_rp, to compute the terminal cost and constraints and the entry points of the OCP for specific values of the decision variables, as well as their first-order and second-order partial derivatives with respect to the decision variables (that is, the gradients and Hessian matrices). The adaptive approaches to deal with plant-model mismatch are directly implemented in each code file where they are used.

Below, the core functions in this folder are described. For more details about their use, check the examples in the files test_hessian_*.m.

The function hessian_compile in the file hessian_compile.m is used to compile the files to be used for numerical integration and computation of the terminal cost and constraints and the entry points of the OCP for specific values of the decision variables, as well as their first-order and second-order partial derivatives with respect to the decision variables, for a specific example. This function [nci,phi,dphidxt,d2phidxt2,f,dfdxp,d2Hdxp2,d2fdxp2,h,dhdx] = hessian_compile(nci,clhs,crhs,lb,ub,phi_form,f_form,gh_form,h_form,np,ni,nt,nti,deriv,varargin) is described as follows:
- nci as an input is a vector that specifies which case is activated when the corresponding state path constraint in h_form becomes active and a switching to a state constraint-seeking arc occurs, and nci as an output also specifies which case is activated when the path constraints specified by lb, ub, and gh_form become active.
- clhs is a cell array with the input names.
- crhs is a cell array where each cell corresponds to a cell array with the values assigned to the inputs in a specific case.
- lb is a vector of lower input bounds.
- ub is a vector of upper input bounds.
- phi_form is a cell array with the formulas of the terminal cost and constraints as a function of the final time t, the final states xn, and the intermediate states xim_n for the mth intermediate time. phi_form may have several columns in the case of terminal cost or constraints that are computed in a loop. In that case, xlk_n denotes the states in the kth loop, ph0_k denotes the cost or constraint in the kth loop, and vk denotes the variable v in the kth loop if those variables are specified in varargin.
- f_form is a cell array with the formulas of the state dynamics as a function of the states xn and the inputs up. f_form may depend on the states xl1_n in a loop if varargin is specified.
- gh_form is a cell array with the formulas of the mixed path constraints.
- h_form is a cell array with the formulas of the state path constraints.
- np is the number of model parameters.
- ni is the number of sensitivity-seeking arcs to be described by additional states.
- nt is the number of switching times including the final time.
- nti is the number of intermediate times used in the computation of the terminal cost and constraints.
- deriv is the higher-order partial derivative of the terminal cost and constraints that needs to be computed.
- If varargin is specified: varargin{1} is a cell array where the first row contains variable names and the second rows contains their values for a loop; varargin{2} is a vector of the length of each loop for computation of the terminal cost and constraints; and varargin{3} is the length of the loop for computation of state dynamics. If varargin is not specified, it is assumed that no loops are used for computing the terminal cost and constraints.
- The outputs phi,dphidxt,d2phidxt2,f,dfdxp,d2Hdxp2,d2fdxp2,h,dhdx are the terminal cost and constraints, state dynamics, and path constraints, as well as their partial derivatives with respect to x, that are computed by the function.

The function hessian_check in the file hessian_check.m is used to compute the terminal cost and constraints and entry points for specified values of the decision variables and to check that various partial derivatives with respect to the decision variables are correctly computed. This function hessian_check(nci,phicase,av,idx_x0,ti,ts,tf,x0,p0,sc,deriv) is described as follows:
- nci is described as the output of hessian_compile with the same name.
- phicase is a vector of indices of the terminal cost and constraints and corresponding partial derivatives that are computed.
- av is a vector that describes the sequence of arcs to be evaluated. For a value less than 0, the arc corresponds to the input constraint-seeking arc specified by the order lb-ub-gh_form of the arguments passed to hessian_compile; for a value greater than 0, the arc corresponds to one of the sensitivity-seeking arcs specified by the argument ni passed to hessian_compile.
- idx_x0 is a vector of indices of x0 for which partial derivatives are computed.
- ti is a vector of intermediate times used in the computation of the terminal cost and constraints.
- ts is a vector of switching times to input constraint-seeking and sensitivity-seeking arcs.
- tf is the final time.
- x0 is a vector of initial states, including the additional states for sensitivity-seeking arcs.
- p0 is a vector of parameters, including the parameters that describe the variation of the additional states for sensitivity-seeking arcs.
- sc is a vector of scaling factors for the decision variables specified by ts, tf, the elements of x0 specified by idx_x0, and p0.
- deriv is described as the input of hessian_compile with the same name.

The function hessian_calc in the file hessian_calc.m is used to compute the terminal cost and constraints, entry points, and other quantities for specified values of the decision variables, as well as the corresponding partial derivatives with respect to those decision variables. This function [J,th,dJdts,dJdtf,dJdx0,dJdp0,dthdts,dthdtf,dthdx0,dthdp0,tlim,xplim,cqlim,drlim,tout,xpout,cqout,drout,uout,d2Jdtsdu,d2Jdtfdu,d2Jdx0du,d2Jdp0du,d2thdtsdu,d2thdtfdu,d2thdx0du,d2thdp0du,dtdulim,dxpdulim,dcqdulim,ddrdulim] = hessian_calc(nci,phicase,av,ti,t0,ts,tf,x0,p0,lam,pi0) is described as follows:
- nci is described as the output of hessian_compile with the same name.
- phicase is described as the input of hessian_check with the same name.
- av is described as the input of hessian_check with the same name.
- ti is described as the input of hessian_check with the same name.
- t0 is the initial time.
- ts is described as the input of hessian_check with the same name.
- tf is described as the input of hessian_check with the same name.
- x0 is described as the input of hessian_check with the same name.
- p0 is described as the input of hessian_check with the same name.
- lam is a matrix of modifiers, where each row contains the modifiers to be added to the corresponding terminal cost or constraint.
- pi0 is a vector of nominal values ts,tf,x0,p0 with respect to which the modifiers are added.
- J is a cell array, where each cell contains the value of a terminal cost or constraint.
- th is a cell array, where each cell contains the value of an entry point in a state constraint-seeking arc.
- dJdts,dJdtf,dJdx0,dJdp0 are cell arrays, where each cell contains the values of the first-order partial derivatives of a terminal cost or constraint with respect to ts,tf,x0,p0, respectively.
- dthdts,dthdtf,dthdx0,dthdp0 are cell arrays, where each cell contains the values of the first-order partial derivatives of an entry point in a state constraint-seeking arc with respect to ts,tf,x0,p0, respectively.
- tlim is a cell array, where the only cell contains a vector with the switching times to all the arcs, the intermediate times, and the final time, in sorted order.
- xplim is a cell array, where the only cell contains a matrix such that each column contains the states and parameters for the corresponding time in tlim.
- cqlim is a cell array, where each cell corresponds to a terminal cost or constraint and contains a matrix such that each column contains the sensitivities with respect to the states and parameters for the corresponding time in tlim.
- drlim is a cell array, where each cell corresponds to an entry point in a state constraint-seeking arc and contains a matrix such that each column contains the sensitivities with respect to the states and parameters for the corresponding time in tlim.
- tout is a vector of times along the trajectory.
- xpout is a matrix such that each column contains the states and parameters for the corresponding time in tlim.
- cqout is a cell array, where each cell corresponds to a terminal cost or constraint and contains a matrix such that each column contains the sensitivities with respect to the states and parameters for the corresponding time in tlim.
- drout is a cell array, where each cell corresponds to an entry point in a state constraint-seeking arc and contains a matrix such that each column contains the sensitivities with respect to the states and parameters for the corresponding time in tlim.
- uout is a matrix such that each column contains the inputs for the corresponding time in tlim.
- d2Jdtsdu,d2Jdtfdu,d2Jdx0du,d2Jdp0du are cell arrays, where each cell contains the second-order partial derivatives of a terminal cost or constraint with respect to ts,tf,x0,p0 and all the decision variables, respectively.
- d2thdtsdu,d2thdtfdu,d2thdx0du,d2thdp0du are cell arrays, where each cell contains the second-order partial derivatives of an entry point in a state constraint-seeking arc with respect to ts,tf,x0,p0 and all the decision variables, respectively.
- dtdulim,dxpdulim,dcqdulim,ddrdulim are cell arrays, where each cell contains the first-order partial derivatives of tlim,xplim,cqlim,drlim with respect to all the decision variables.

The function hessian_solve in the file hessian_solve.m is used to solve a single-input OCP to global optimality. The function [tau,p_sol,d_sol,deg,tel,fvalv,uv,N,p_mon,p_s,d,dus,flag,cphi,As,ys,y,Av,yv,u0,du,sc] = hessian_solve(tf,x0,nci,nphi,av,np,ni,nt,ti,scv,deriv,lb,ub,varargin) is described as follows:
- tf is described as the input of hessian_check with the same name.
- x0 is a vector of initial states, excluding the additional states for sensitivity-seeking arcs.
- nci is described as the output of hessian_compile with the same name.
- nphi is the number of terminal cost and constraints.
- av is either described as the input of hessian_check with the same name, or the number of arcs in the arc sequences to be investigated.
- np is described as the input of hessian_compile with the same name.
- ni is described as the input of hessian_compile with the same name.
- nt is the number of switching times including the final time only if it is a decision variable.
- ti is described as the input of hessian_check with the same name.
- scv is a vector of scaling factors to be used by hessian_fmincon for the decision variables that correspond to the nt switching times, the ni additional states, and ni additional parameters if deriv is different from 0.
- deriv is the scaling factor for the ni additional parameters if deriv is different from 0, or it indicates that the ni additional parameters should be computed from an initial local solution if deriv is equal to 0.
- lb is described as the input of hessian_compile with the same name.
- ub is described as the input of hessian_compile with the same name.
- If varargin is specified: varargin{1} is the degree of the polynomials used to approximate the terminal cost and constraints; varargin{2} is the number of points used for this polynomial approximation; varargin{3} specifies whether Chebyshev nodes are used to choose the points for polynomial approximation; varargin{4} specifies whether supplementary points are computed in a sub-region of the space of decision variables; varargin{5} is described as the input lam of hessian_calc; and varargin{6} is described as the input pi0 of hessian_calc. If varargin is not specified, default values are used.
- with the exception of u0,du,sc, all the outputs are cell arrays, where each cell contains outputs for a specific arc sequence.
- tau{l} is the globally optimal cost of the POP computed from the solution to the SDP.
- p_sol{l} is the globally optimal solution to the POP computed from the primal solution to the SDP.
- d_sol{l} is the globally optimal solution to the POP computed from the dual solution to the SDP.
- deg{l} is the relaxation order of the SDP for which a globally optimal solution to the POP is computed.
- tel{l} is the elapsed time.
- fvalv{l} is the globally optimal cost of the OCP for the arc sequence computed from the solution to the SDP followed by local optimization of the OCP without approximation via polynomials.
- uv{l} is the globally optimal solution to the OCP for the arc sequence computed from the solution to the SDP followed by local optimization of the OCP without approximation via polynomials.
- N{l} is the number of decision variables.
- p_mon{l} is a matrix, where each row contains the powers of a monomial used for the polynomial approximation.
- p_s{l} is the number of monomials used for the polynomial approximation.
- d{l} is a matrix, where each column contains the polynomial coefficients of an inequality constraint that avoids inadmissible values of decision variables.
- dus{l} is a matrix, where each row contains a sample point in the space of decision variables in terms of offsets with respect to their nominal values.
- flag{l} is a vector that specifies which sample points in the rows of dus{l} are admissible for computation of polynomial coefficients.
- cphi{l} is a cell array, where each cell contains the polynomial coefficients for a terminal cost or constraint and partial derivatives.
- As{l} is a matrix, where each row contains the scaled regressors for an admissible sample point in the space of decision variables.
- ys{l} is a cell array, where each cell contains the values of a terminal cost or constraint and partial derivatives for the sample points in the rows of dus{l}.
- y{l} is a cell array, where each cell contains the values of a terminal cost or constraint and partial derivatives for the globally optimal solution in uv{l}.
- Av{l} is a row vector that contains the scaled regressors for the globally optimal solution in uv{l}.
- yv{l} is a cell array, where each cell contains a vector with the values of the polynomial approximations of the terminal cost or constraint and partial derivatives for the globally optimal solution in uv{l}.
- u0 is a vector with the nominal values of the decision variables.
- du is a two-row matrix, where each column specifies the maximum negative and positive offsets used for each decision variable.
- sc is a vector with the scaling factors used for the decision variables.

The function hessian_fmincon in the file hessian_fmincon.m is used to solve a single-input OCP to local optimality from an initial guess via MATLAB's fmincon. The function [u,tau,exitflag,output,nu,signth] = hessian_fmincon(nci,nphi,av,np,idx_x0,ti,ts0,tf0,x00,p00,sc,lam,pi0,lb,ub,ibp) is described as follows:
- nci is described as the output of hessian_compile with the same name.
- nphi is described as the input of hessian_solve with the same name.
- av is described as the input of hessian_check with the same name.
- np is described as the input of hessian_compile with the same name.
- idx_x0 is a vector of indices of x00 that are considered as decision variables.
- ti is described as the input of hessian_check with the same name.
- ts0 is the initial guess for the vector of switching times to input constraint-seeking and sensitivity-seeking arcs.
- tf is the initial guess for the final time.
- x0 is the initial guess for the vector of initial states, including the additional states for sensitivity-seeking arcs.
- p0 is the initial guess for the vector of parameters, including the parameters that describe the variation of the additional states for sensitivity-seeking arcs.
- sc is a vector of scaling factors for the decision variables specified by ts0, tf0, the elements of x00 specified by idx_x0, and p00.
- lam is described as the input of hessian_calc with the same name.
- pi0 is described as the input of hessian_calc with the same name.
- lb is described as the input of hessian_compile with the same name.
- ub is described as the input of hessian_compile with the same name.
- ibp is a parameter that describes how aggressively a local solution is searched from the initial guess.
- u is a locally optimal solution to the OCP.
- tau is a locally optimal cost of the OCP.
- exitflag,output,nu are similar to the outputs of MATLAB's fmincon with the same names.
- signth is a vector that specifies the sign of the last entry point with respect to the optimal switching times and final time.

The function hessian_plot in the file hessian_plot.m is used to plot the states, inputs, and sensitivities with respect to the states and parameters for specific values of the decision variables. The function hessian_plot(nci,phicase,av,ti,t0,ts,tf,x0,p0,lam,pi0,nu_ineq) is described by the same input arguments as hessian_calc, and in addition nu_ineq specifies the Lagrange multipliers for the terminal constraints.
