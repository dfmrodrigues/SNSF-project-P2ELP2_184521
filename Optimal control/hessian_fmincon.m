function [u,tau,exitflag,output,nu,signth] = hessian_fmincon(nci,nphi,av,np,idx_x0,ti,ts0,tf0,x00,p00,sc,lam,pi0,lb,ub,ibp)
phicase = 1:nphi;
nt = length(ts0)+1;
nix = length(idx_x0);
np = nix+np;
last_du = [];
last_obj = [];
last_ineq = [];
last_gobj = [];
last_gineq = [];
ss = av>0;
du0 = zeros(1,nt+nix+np);
A = [[zeros(1,nt);eye(nt-1,nt)/sc(nt);zeros(4*nix,nt)]-[eye(nt)/sc(nt);zeros(4*nix,nt)],...
    [zeros(nt,2*nix);eye(2*nix)./[sc(nt+(1:nix))';sc(nt+(1:nix))'];-eye(2*nix)./[sc(nt+(1:nix))';sc(nt+(1:nix))']]];
b = [zeros(nt,1);ub*ones(2*nix,1)./[sc(nt+(1:nix))';sc(nt+(1:nix))'];-lb*ones(2*nix,1)./[sc(nt+(1:nix))';sc(nt+(1:nix))']];
options = optimoptions('fmincon','Algorithm','interior-point','Display','iter-detailed',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',false,...
    'FunctionTolerance',1e-12,'OptimalityTolerance',1e-12,'ConstraintTolerance',1e-9,'StepTolerance',1e-15,...
    'MaxFunctionEvaluations',3000,'InitBarrierParam',ibp);
addpath('./hessian');
[du,tau,exitflag,output,nu] = fmincon(@objective,du0,[],[],[],[],[],[],@constraints,options);
u = [ts0,tf0,x00(idx_x0)',p00]+sc.*du;
rmpath('./hessian');
    function change2new(du)
        last_du = du;
        ts = ts0+sc(1:(nt-1)).*du(1:(nt-1));
        tf = tf0+sc(nt).*du(nt);
        x0 = x00;
        x0(idx_x0) = x0(idx_x0)+(sc(nt+(1:nix)).*du(nt+(1:nix)))';
        p0 = p00'+(sc(nt+nix+(1:np)).*du(nt+nix+(1:np)))';
        tb = [0,ts];
        te = [ts,tf];
        xf = x0(idx_x0)+diag(p0)*(te(ss)-tb(ss))';
        dxfdte = zeros(np,nt);
        dxfdte(:,ss) = diag(p0);
        dxfdtb = zeros(np,nt);
        dxfdtb(:,ss) = -diag(p0);
        if(prod(A(1:nt,:)*[ts,tf,x0(idx_x0)',xf']'-b(1:nt)<=-1e-9)==0||...
                (~isempty(last_obj)&&prod(A((nt+1):end,:)*[ts,tf,x0(idx_x0)',xf']'-b((nt+1):end)<=-1e-9)==0))
            last_obj = Inf;
            last_gobj = zeros(nt+nix+np,1);
            last_ineq = zeros(nphi-1+nt+4*nix,1);
            last_gineq = zeros(nt+nix+np,nphi-1+nt+4*nix);
        else
            [J,th,dJdts,dJdtf,dJdx0,dJdp0] = hessian_calc(nci,phicase,av,ti,0,ts,tf,x0,p0,lam,pi0);
            if(isempty(th))
                signth = ones(1,nt);
            else
                signth = sign(cell2mat(th(length(th)))-[ts,tf]);
            end
            obj = cell2mat(J(1)');
            dobjdts = cell2mat(dJdts(1)');
            dobjdtf = cell2mat(dJdtf(1)');
            dobjdx0 = cell2mat(dJdx0(1)');
            dobjdp0 = cell2mat(dJdp0(1)');
            dobjdx0 = dobjdx0(:,idx_x0);
            if(nphi==1)
                ineq = zeros(0,size(obj,2));
                dineqdts = zeros(0,size(dobjdts,2));
                dineqdtf = zeros(0,size(dobjdtf,2));
                dineqdx0 = zeros(0,size(dobjdx0,2));
                dineqdp0 = zeros(0,size(dobjdp0,2));
            else
                ineq = cell2mat(J(2:nphi)');
                dineqdts = cell2mat(dJdts(2:nphi)');
                dineqdtf = cell2mat(dJdtf(2:nphi)');
                dineqdx0 = cell2mat(dJdx0(2:nphi)');
                dineqdp0 = cell2mat(dJdp0(2:nphi)');
                dineqdx0 = dineqdx0(:,idx_x0);
            end
            sens = 1;
            ineq = [ineq;sens*(A*[ts,tf,x0(idx_x0)',xf']'-b)];
            dineqdts = [dineqdts;sens*(A(:,1:(nt-1))+A(:,nt+nix+(1:np))*(dxfdtb(:,2:end)+dxfdte(:,1:end-1)))];
            dineqdtf = [dineqdtf;sens*(A(:,nt)+A(:,nt+nix+(1:np))*dxfdte(:,end))];
            dineqdx0 = [dineqdx0;sens*(A(:,nt+(1:nix))+A(:,nt+nix+(1:np)))];
            dineqdp0 = [dineqdp0;sens*A(:,nt+nix+(1:np)).*(te(ss)-tb(ss))];
            last_obj = obj;
            last_ineq = ineq;
            last_gobj = [sc(1:(nt-1)).*dobjdts,sc(nt).*dobjdtf,sc(nt+(1:nix)).*dobjdx0,sc(nt+nix+(1:np)).*dobjdp0]';
            last_gineq = [sc(1:(nt-1)).*dineqdts,sc(nt).*dineqdtf,sc(nt+(1:nix)).*dineqdx0,sc(nt+nix+(1:np)).*dineqdp0]';
        end
    end
    function [obj,gobj] = objective(du)
        if(~isequal(du,last_du))
            change2new(du);
        end
        obj = last_obj;
        gobj = last_gobj;
    end
    function [ineq,eq,gineq,geq] = constraints(du)
        if(~isequal(du,last_du))
            change2new(du);
        end
        ineq = last_ineq;
        gineq = last_gineq;
        eq = [];
        geq = [];
    end
end