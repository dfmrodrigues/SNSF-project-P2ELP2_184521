function [wopt,thetaopt,pols,dpols,innerp,diff_pols_eval,diff_pols_deval,diff_basis_eval,diff_basis_deval] = quadrature(nparam,max_order,orth_pols,pdf,w,theta)
ords = calc_mon(max_order,nparam);
form = cell(max_order+1,nparam);
dform = cell(max_order+1,nparam);
pols = cell(max_order+1,nparam);
dpols = cell(max_order+1,nparam);
innerp = zeros(max_order+1,max_order+1,nparam);
fidp = fopen('pols.c','w');
fprintf(fidp,'#include <mex.h>\n#include <math.h>\n\n');
fprintf(fidp,'void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])\n{\n');
fprintf(fidp,'double *pol, *dpol, *x;\n\n');
fprintf(fidp,'x = mxGetPr(prhs[0]);\n\n');
fprintf(fidp,'plhs[0] = mxCreateDoubleMatrix(%d,%d,mxREAL);\n',nparam,max_order+1);
fprintf(fidp,'pol = mxGetPr(plhs[0]);\n\n');
fprintf(fidp,'plhs[1] = mxCreateDoubleMatrix(%d,%d,mxREAL);\n',nparam,max_order+1);
fprintf(fidp,'dpol = mxGetPr(plhs[1]);\n\n');
for k = 1:nparam
    for n = 0:max_order
        form{n+1,k} = orth_pols(k,n);
        dform{n+1,k} = diff(form{n+1,k},sym('x'));
        formchar = ccode(form{n+1,k});
        formchar = replace(formchar,'  t0 = ','');
        formchar = replace(formchar,'x',['x[',num2str(k-1),']']);
        dformchar = ccode(dform{n+1,k});
        dformchar = replace(dformchar,'  t0 = ','');
        dformchar = replace(dformchar,'x',['x[',num2str(k-1),']']);
        fprintf(fidp,['pol[%d] = ',formchar,'\n'],k-1+nparam*n);
        fprintf(fidp,['dpol[%d] = ',dformchar,'\n'],k-1+nparam*n);
        pols{n+1,k} = str2func(['@(x)',replace(char(form{n+1,k}),'^','.^')]);
        dpols{n+1,k} = str2func(['@(x)',replace(char(dform{n+1,k}),'^','.^')]);
    end
    for n1 = 0:max_order
        for n2 = 0:max_order
            innerp(n1+1,n2+1,k) = integral(@(x)pols{n1+1,k}(x).*pols{n2+1,k}(x).*pdf(x,k),-Inf,Inf,'AbsTol',1e-7);
        end
    end
end
fprintf(fidp,'}\n');
fclose(fidp);
eval('mex -output pols_mex pols.c');
eval('mex -output basis_mex basis.c');

mtrain = size(theta,1);
points = reshape(theta,[mtrain*nparam,1]);
if(nparam==2)
    [X,Y] = meshgrid(-1:0.01:1,-1:0.01:1);
    Z = pdf(X,1).*pdf(Y,2);
    figure,contourf(X,Y,Z,100,'LineStyle','none'),colormap jet,set(gca,'XLim',[-1,1],'YLim',[-1,1]),hold on;
end
x = [w;points];

options = optimoptions('fmincon','Algorithm','trust-region-reflective','Display','iter-detailed','MaxFunctionEvaluations',20000,'MaxIterations',20000,...
    'OptimalityTolerance',eps^2,'FunctionTolerance',eps^2,'StepTolerance',eps,'SpecifyObjectiveGradient',true);
xopt = fmincon(@(x)quad_obj(x,mtrain,nparam,max_order,ords,pols,dpols,true),x,[],[],[],[],[zeros(mtrain,1);-Inf*ones(mtrain*nparam,1)],[],[],options);
wopt = xopt(1:mtrain);
pointsopt = xopt(mtrain+1:end);
% options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter','MaxFunctionEvaluations',20000,'MaxIterations',20000,...
%     'OptimalityTolerance',eps,'StepTolerance',eps,'SpecifyObjectiveGradient',true,'CheckGradient',true);
% xopt = lsqnonlin(@(x)quad_obj(x,mtrain,nparam,max_order,ords,pols,dpols,false),x,[],[],options);
% wopt = xopt(1:mtrain);
% pointsopt = xopt(mtrain+1:end);
% [~,~,C,d] = quad_obj(x,mtrain,nparam,max_order,ords,pols,dpols,false);
% [xopt,resopt] = lsqnonneg(C,d);
% resopt
% wopt = xopt;
% pointsopt = points;
thetaopt = reshape(pointsopt,mtrain,nparam);
if(nparam==2)
    scatter(thetaopt(:,1),thetaopt(:,2),max(wopt,1e-6)*1000,'ko','LineWidth',2),hold on;
    scatter(thetaopt(:,1),thetaopt(:,2),max(wopt,1e-6)*1000,'wx','LineWidth',2);
    set(gca,'TickLabelInterpreter','latex','FontSize',15);
    hXL = xlabel('$\theta_1$','Interpreter','latex');
    set(hXL,'FontSize',25);
    hYL = ylabel('$\theta_2$','Interpreter','latex');
    set(hYL,'FontSize',25);
end
[~,~,~,~,diff_pols_eval,diff_pols_deval,diff_basis_eval,diff_basis_deval] = quad_obj(x,mtrain,nparam,max_order,ords,pols,dpols,false);

save('quadrature.mat','wopt','thetaopt');
end