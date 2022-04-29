function [nci,phi,dphidxt,d2phidxt2,f,dfdxp,d2Hdxp2,d2fdxp2,h,dhdx] = hessian_compile(nci,clhs,crhs,lb,ub,phi_form,f_form,gh_form,h_form,np,ni,nt,nti,deriv,varargin)
if(isempty(varargin))
    vars = {};
    lims_phi = [];
    lims_f = 1;
else
    vars = varargin{1};
    lims_phi = varargin{2};
    lims_f = varargin{3};
end
nf = length(f_form);
nc = length(crhs);
nu = length(clhs);
nci = [nci+ni,-(1:(length(lb)+length(ub)+length(gh_form)))];
ci = [1:ni,ni+(1:nc),-(1:(length(lb)+length(ub)+length(gh_form)))];
crhs = [cell(1,ni),crhs,cell(1,length(lb)+length(ub)+length(gh_form))];
prhs = cell(1,ni+nc+length(lb)+length(ub)+length(gh_form));
for i = 1:ni
    s = cell(nu,1);
    for k = 1:nu
        s(k) = {['x',num2str(nf*lims_f+(i-1)*nu+k)]};
    end
    crhs(:,i) = {str2sym(s)};
    s = cell(nu*ni,1);
    for j = 1:ni
        for k = 1:nu
            if(i==j)
                s((j-1)*nu+k) = {['p',num2str(np+(i-1)*nu+k)]};
            else
                s((j-1)*nu+k) = {'0'};
            end
        end
    end
    prhs(:,i) = {str2sym(s)};
end
for i = 1:nc
    s = cell(nu*ni,1);
    for j = 1:ni
        for k = 1:nu
            s((j-1)*nu+k) = {'0'};
        end
    end
    prhs(:,ni+i) = {str2sym(s)};
end
for i = 1:length(lb)
    s = cell(nu,1);
    for j = 1:nu
        if(i==j)
            s(j) = {num2str(lb(i))};
        else
            s(j) = {'0'};
        end
    end
    crhs(:,ni+nc+i) = {str2sym(s)};
    s = cell(nu*ni,1);
    for j = 1:ni
        for k = 1:nu
            s((j-1)*nu+k) = {'0'};
        end
    end
    prhs(:,ni+nc+i) = {str2sym(s)};
end
s = cell(length(lb),1);
for i = 1:length(lb)
    s(i) = {str2sym([num2str(lb(i)),'-',char(clhs(i)),'-1E-8'])};
end
h_form = [h_form;s];
for i = 1:length(ub)
    s = cell(nu,1);
    for j = 1:nu
        if(i==j)
            s(j) = {num2str(ub(i))};
        else
            s(j) = {'0'};
        end
    end
    crhs(:,ni+nc+length(lb)+i) = {str2sym(s)};
    s = cell(nu*ni,1);
    for j = 1:ni
        for k = 1:nu
            s((j-1)*nu+k) = {'0'};
        end
    end
    prhs(:,ni+nc+length(lb)+i) = {str2sym(s)};
end
s = cell(length(ub),1);
for i = 1:length(ub)
    s(i) = {str2sym([char(clhs(i)),'-',num2str(ub(i)),'-1E-8'])};
end
h_form = [h_form;s];
for i = 1:length(gh_form)
    ghfun = gh_form(i);
    var = [];
    for l = 0:1
        s = cell(nu,1);
        for j = 1:nu
            if(isempty(var))
                var = solve(ghfun,clhs(j));
                if(~isempty(var))
                    s(j) = {char(var)};
                    continue;
                end
            end
            s(j) = {'0'};
        end
        if(~isempty(var))
            break;
        else
            ghfun = jacobian(ghfun,sym('x',[1,length(f_form)]))*f_form;
        end
    end
    crhs(:,ni+nc+length(lb)+length(ub)+i) = {str2sym(s)};
    s = cell(nu*ni,1);
    for j = 1:ni
        for k = 1:nu
            s((j-1)*nu+k) = {'0'};
        end
    end
    prhs(:,ni+nc+length(lb)+length(ub)+i) = {str2sym(s)};
end
s = cell(length(gh_form),1);
for i = 1:length(gh_form)
    s(i) = {gh_form(i)};
end
h_form = [h_form;s];
nh = length(h_form);
nx = nf*lims_f+ni*nu;
np = np+ni*nu;
deg = 1;
u_mon = calc_mon(deg,nx+nt+np);
u_s = size(u_mon,1);
t = sym('t');
if(deriv>=2)
    x = sym('x',[1,nx*u_s]);
    c = sym('c',[1,nx*u_s]);
    xi = sym('xi',[nti,nx*u_s]);
    xl = sym('xl%d_%d',[length(lims_phi),nf*u_s]);
else
    x = sym('x',[1,nx]);
    c = sym('c',[1,nx]);
    xi = sym('xi',[nti,nx]);
    xl = sym('xl%d_%d',[length(lims_phi),nf]);
end
p = sym('p',[1,np]);
phi = cell(size(phi_form));
phichar = cell(size(phi_form));
f = cell(size(crhs));
fchar = cell(size(crhs));
u = cell(size(crhs));
uchar = cell(size(crhs));
if(deriv>=1)
    dphidxt = cell(size(phi_form));
    dphidxtchar = cell(size(phi_form));
    dfdxp = cell(size(crhs));
    dfdxpchar = cell(size(crhs));
    dudxp = cell(size(crhs));
    dudxpchar = cell(size(crhs));
end
if(deriv>=2)
    d2phidxt2 = cell(size(phi_form));
    d2phidxt2char = cell(size(phi_form));
    d2Hdxp2 = cell(size(crhs));
    d2fdxp2 = cell(size(crhs));
    d2fdxp2char = cell(size(crhs));
    dndydtdunchar = cell(size(crhs));
    dndadtdunchar = cell(size(crhs));
end
delete('hessian/*_hessian.c');
delete('hessian/*_hessian_mex.mexmaci64');
fidp0 = fopen('hessian/phixt0_hessian.c','w');
fidy0 = fopen('hessian/f0_hessian.c','w');
fida0 = fopen('hessian/f0x_hessian.c','w');
if(deriv>=1)
    fidp1 = fopen('hessian/phixt1_hessian.c','w');
end
if(deriv>=2)
    fidp2 = fopen('hessian/phixt2_hessian.c','w');
    fidy1 = fopen('hessian/f1_hessian.c','w');
    fida1 = fopen('hessian/f1x_hessian.c','w');
end
fidgy = fopen('hessian/gh_hessian.c','w');
fidga = fopen('hessian/ghx_hessian.c','w');
fprintf(fidp0,'#include <mex.h>\n#include <math.h>\n\n');
fprintf(fidp0,'void phixt0fun(double *phi_xt0, double t, double *x, double *xi, const mxArray *args[])\n{\n');
fprintf(fidp0,'int cases;\n\n');
for l = 1:size(vars,2)
    fprintf(fidp0,['const double ',vars{1,l},'[',num2str(length(vars{2,l})),'] = {',strip(replace(sprintf('%15.15g\n',vars{2,l}{:}),newline,','),'right',','),'};\n\n']);
end
for k = 2:size(phi_form,2)
    fprintf(fidp0,['int l',num2str(k-1),';\n\n']);
end
for k = 2:size(phi_form,2)
    fprintf(fidp0,['double phi_xt0_',num2str(k-1),'[1];\n\n']);
end
fprintf(fidp0,'cases = mxGetScalar(args[0]);\n\n');
fprintf(fidp0,'switch(cases) {\n');
fprintf(fidy0,'#include <mex.h>\n#include <math.h>\n\n');
fprintf(fidy0,'void ffun(double *f, double *u, double *x, const mxArray *args[])\n{\n');
fprintf(fidy0,'int cases;\ndouble *p = x+%d;\n\n',nx);
for l = 1:size(vars,2)
    fprintf(fidy0,['const double ',vars{1,l},'[',num2str(length(vars{2,l})),'] = {',strip(replace(sprintf('%15.15g\n',vars{2,l}{:}),newline,','),'right',','),'};\n\n']);
end
fprintf(fidy0,'int l1;\n\n');
fprintf(fidy0,'cases = mxGetScalar(args[0]);\n\n');
fprintf(fidy0,'switch(cases) {\n');
fprintf(fida0,'#include <mex.h>\n#include <math.h>\n\n');
fprintf(fida0,'void fxfun(double *f_x, double *u_x, double *x, const mxArray *args[])\n{\n');
fprintf(fida0,'int cases;\ndouble *p = x+%d;\n\n',nx);
if(deriv>=1)
    fprintf(fidp1,'#include <mex.h>\n#include <math.h>\n\n');
    fprintf(fidp1,'void phixt1fun(double *phi_xt1, double t, double *x, double *xi, const mxArray *args[])\n{\n');
    fprintf(fidp1,'int cases;\n\n');
    for l = 1:size(vars,2)
        fprintf(fidp1,['const double ',vars{1,l},'[',num2str(length(vars{2,l})),'] = {',strip(replace(sprintf('%15.15g\n',vars{2,l}{:}),newline,','),'right',','),'};\n\n']);
    end
    for k = 2:size(phi_form,2)
        fprintf(fidp1,['int l',num2str(k-1),';\n\n']);
    end
    for k = 2:size(phi_form,2)
        fprintf(fidp1,['double phi_xt0_',num2str(k-1),'[1];\n\n']);
        fprintf(fidp1,['double phi_xt1_',num2str(k-1),'[',num2str(nx*(1+nti)+1),'];\n\n']);
    end
    fprintf(fidp1,'cases = mxGetScalar(args[0]);\n\n');
    fprintf(fidp1,'switch(cases) {\n');
    for l = 1:size(vars,2)
        fprintf(fida0,['const double ',vars{1,l},'[',num2str(length(vars{2,l})),'] = {',strip(replace(sprintf('%15.15g\n',vars{2,l}{:}),newline,','),'right',','),'};\n\n']);
    end
    fprintf(fida0,'int l1;\n\n');
end
fprintf(fida0,'cases = mxGetScalar(args[0]);\n\n');
fprintf(fida0,'switch(cases) {\n');
if(deriv>=2)
    fprintf(fidp2,'#include <mex.h>\n#include <math.h>\n\n');
    fprintf(fidp2,'void phixt2fun(double *phi_xt2, double t, double *x, double *xi, const mxArray *args[])\n{\n');
    fprintf(fidp2,'int cases;\n\n');
    for l = 1:size(vars,2)
        fprintf(fidp2,['const double ',vars{1,l},'[',num2str(length(vars{2,l})),'] = {',strip(replace(sprintf('%15.15g\n',vars{2,l}{:}),newline,','),'right',','),'};\n\n']);
    end
    for k = 2:size(phi_form,2)
        fprintf(fidp2,['int l',num2str(k-1),';\n\n']);
    end
    for k = 2:size(phi_form,2)
        fprintf(fidp2,['double phi_xt0_',num2str(k-1),'[1];\n\n']);
        fprintf(fidp2,['double phi_xt2_',num2str(k-1),'[',num2str((nx*(1+nti)+1)*(nx*(1+nti)+1)),'];\n\n']);
    end
    fprintf(fidp2,'cases = mxGetScalar(args[0]);\n\n');
    fprintf(fidp2,'switch(cases) {\n');
    fprintf(fidy1,'#include <mex.h>\n#include <math.h>\n\n');
    fprintf(fidy1,'void ffun(double *f, double *u, double *x, const mxArray *args[])\n{\n');
    fprintf(fidy1,'int cases;\ndouble *p = x+%d;\n\n',nx);
    fprintf(fidy1,'int l1;\n\n');
    fprintf(fidy1,'cases = mxGetScalar(args[0]);\n\n');
    fprintf(fidy1,'switch(cases) {\n');
    fprintf(fida1,'#include <mex.h>\n#include <math.h>\n\n');
    fprintf(fida1,'void fxfun(double *f_x, double *u_x, double *x, const mxArray *args[])\n{\n');
    fprintf(fida1,'int cases;\ndouble *p = x+%d;\n\n',nx);
    fprintf(fida1,'int l1;\n\n');
    fprintf(fida1,'cases = mxGetScalar(args[0]);\n\n');
    fprintf(fida1,'switch(cases) {\n');
end
fprintf(fidgy,'#include <mex.h>\n#include <math.h>\n\n');
fprintf(fidgy,'void ghfun(double *gh, double *u, double *x, const mxArray *args[])\n{\n');
fprintf(fidga,'#include <mex.h>\n#include <math.h>\n\n');
fprintf(fidga,'void ghxfun(double *gh_x, double *u_x, double *x, const mxArray *args[])\n{\n');
for i = 1:size(phi_form,1)
for k = 1:size(phi_form,2)
    phi{i,k} = phi_form(i,k);
    phichar{i,k} = arrayfun(@ccode,formula(phi{i,k}),'Uniform',0);
    phichar{i,k} = replace(phichar{i,k},'  t0 = ','');
    phichar{i,k} = replace(phichar{i,k},';','');
    for j = 1:(k-1)
        for l = nf:-1:1
            phichar{i,k} = replace(phichar{i,k},['x','l',num2str(j),'_',num2str(l)],['x[',num2str(l-1),'+(l',num2str(j),'-1)*',num2str(nf),']']);
        end
    end
    for l = nx:-1:1
        phichar{i,k} = replace(phichar{i,k},['x',num2str(l)],['x[',num2str(l-1),']']);
    end
    for l = (nx*nti):-1:1
        phichar{i,k} = replace(phichar{i,k},['xi',num2str(floor((l-1)/nx)+1),'_',num2str(rem(l-1,nx)+1)],['xi[',num2str(l-1),']']);
    end
    for j = 1:(k-1)
        for l = 1:size(vars,2)
            phichar{i,k} = replace(phichar{i,k},[vars{1,l},'l',num2str(j)],[vars{1,l},'[l',num2str(j),'-1]']);
        end
    end
    phichar{i,k} = replace(phichar{i,k},['ph0_',num2str(k)],['phi_xt0_',num2str(k),'[',num2str(0),']']);
    if(~strcmp(phichar{i,k}{1},'0.0'))
        if(k>1)
            phichar{i,k}{1} = ['phi_xt0_',num2str(k-1),'[',num2str(0),'] += ',phichar{i,k}{1},';'];
        else
            phichar{i,k}{1} = ['phi_xt0[',num2str(0),'] = ',phichar{i,k}{1},';'];
        end
    else
        phichar{i,k}{1} = [];
    end
    if(deriv>=1)
        if(k==size(phi,2)||strcmp(char(formula(phi_form(i,k+1))),'0'))
            dphidxt{i,k} = jacobian(phi{i,k},[t,x(1:nx),reshape(xi(:,1:nx).',[],1).',reshape(xl(:,1:nf).',[],1).']);
        else
            dphidxt{i,k} = repmat(jacobian(phi{i,k},sym(['ph0_',num2str(k)]))*sym(['ph1_',num2str(k)]),1,nx*(1+nti)+1+nf*length(lims_phi));
        end
        dphidxtchar{i,k} = arrayfun(@ccode,formula(dphidxt{i,k}),'Uniform',0);
        dphidxtchar{i,k} = replace(dphidxtchar{i,k},'  t0 = ','');
        dphidxtchar{i,k} = replace(dphidxtchar{i,k},';','');
        for j = 1:(k-1)
            for l = nf:-1:1
                dphidxtchar{i,k} = replace(dphidxtchar{i,k},['x','l',num2str(j),'_',num2str(l)],['x[',num2str(l-1),'+(l',num2str(j),'-1)*',num2str(nf),']']);
            end
        end
        for l = nx:-1:1
            dphidxtchar{i,k} = replace(dphidxtchar{i,k},['x',num2str(l)],['x[',num2str(l-1),']']);
        end
        for l = (nx*nti):-1:1
            dphidxtchar{i,k} = replace(dphidxtchar{i,k},['xi',num2str(floor((l-1)/nx)+1),'_',num2str(rem(l-1,nx)+1)],['xi[',num2str(l-1),']']);
        end
        for j = 1:(k-1)
            for l = 1:size(vars,2)
                dphidxtchar{i,k} = replace(dphidxtchar{i,k},[vars{1,l},'l',num2str(j)],[vars{1,l},'[l',num2str(j),'-1]']);
            end
        end
        dphidxtchar{i,k} = replace(dphidxtchar{i,k},['ph0_',num2str(k)],['phi_xt0_',num2str(k),'[',num2str(0),']']);
        for j = 1:(nx*(1+nti)+1)
            dphidxtchar{i,k}{j} = replace(dphidxtchar{i,k}{j},['ph1_',num2str(k)],...
                ['phi_xt1_',num2str(k),'[',num2str(j-1),']']);
            if(~strcmp(dphidxtchar{i,k}{j},'0.0'))
                if(k>1)
                    dphidxtchar{i,k}{j} = ['phi_xt1_',num2str(k-1),'[',num2str(j-1),'] += ',...
                        dphidxtchar{i,k}{j},';'];
                else
                    dphidxtchar{i,k}{j} = ['phi_xt1[',num2str(j-1),'] = ',...
                        dphidxtchar{i,k}{j},';'];
                end
            else
                dphidxtchar{i,k}{j} = [];
            end
        end
        for j = (nx*(1+nti)+2):(nx*(1+nti)+1+nf*length(lims_phi))
            jl = j-(nx*(1+nti)+1);
            if(~strcmp(dphidxtchar{i,k}{j},'0.0')&&k==size(phi,2))
                dphidxtchar{i,k}{j} = ['phi_xt1_',num2str(k-1),'[',num2str(rem(jl-1,nf)+1),'+(l',num2str(floor((jl-1)/nf)+1),'-1)*',num2str(nf),'] += ',...
                    dphidxtchar{i,k}{j},';'];
            else
                dphidxtchar{i,k}{j} = [];
            end
        end
    else
        dphidxt{i,k} = [];
    end
    if(deriv>=2)
        if(k==size(phi,2)||strcmp(char(formula(phi_form(i,k+1))),'0'))
            d2phidxt2{i,k} = jacobian(dphidxt{i,k}.',[t,x(1:nx),reshape(xi(:,1:nx).',[],1).']);
        else
            d2phidxt2{i,k} = repmat(jacobian(phi{i,k},sym(['ph0_',num2str(k)]))*sym(['ph2_',num2str(k)]),nx*(1+nti)+1,nx*(1+nti)+1);
        end
        d2phidxt2char{i,k} = arrayfun(@ccode,formula(d2phidxt2{i,k}),'Uniform',0);
        d2phidxt2char{i,k} = replace(d2phidxt2char{i,k},'  t0 = ','');
        d2phidxt2char{i,k} = replace(d2phidxt2char{i,k},';','');
        for l = nx:-1:1
            d2phidxt2char{i,k} = replace(d2phidxt2char{i,k},['x',num2str(l)],['x[',num2str(l-1),']']);
        end
        for l = (nx*nti):-1:1
            d2phidxt2char{i,k} = replace(d2phidxt2char{i,k},['xi',num2str(floor((l-1)/nx)+1),'_',num2str(rem(l-1,nx)+1)],['xi[',num2str(l-1),']']);
        end
        for j = 1:(k-1)
            for l = 1:size(vars,2)
                d2phidxt2char{i,k} = replace(d2phidxt2char{i,k},[vars{1,l},'l',num2str(j)],[vars{1,l},'[l',num2str(j),'-1]']);
            end
        end
        d2phidxt2char{i,k} = replace(d2phidxt2char{i,k},['ph0_',num2str(k)],['phi_xt0_',num2str(k),'[',num2str(0),']']);
        for j = 1:(nx*(1+nti)+1)
            for l = 1:(nx*(1+nti)+1)
                d2phidxt2char{i,k}{j,l} = replace(d2phidxt2char{i,k}{j,l},['ph2_',num2str(k)],...
                    ['phi_xt2_',num2str(k),'[',num2str(j-1+(l-1)*(nx*(1+nti)+1)),']']);
                if(~strcmp(d2phidxt2char{i,k}{j,l},'0.0'))
                    if(k>1)
                        d2phidxt2char{i,k}{j,l} = ['phi_xt2_',num2str(k-1),'[',num2str(j-1+(l-1)*(nx*(1+nti)+1)),'] += ',...
                            d2phidxt2char{i,k}{j,l},';'];
                    else
                        d2phidxt2char{i,k}{j,l} = ['phi_xt2[',num2str(j-1+(l-1)*(nx*(1+nti)+1)),'] = ',...
                            d2phidxt2char{i,k}{j,l},';'];
                    end
                else
                    d2phidxt2char{i,k}{j,l} = [];
                end
            end
        end
    else
        d2phidxt2{i,k} = [];
    end
end
fprintf(fidp0,'case %d:\n',i);
if(deriv>=1)
    fprintf(fidp1,'case %d:\n',i);
end
if(deriv>=2)
    fprintf(fidp2,'case %d:\n',i);
end
for k = 2:size(phi_form,2)
    fprintf(fidp0,['phi_xt0_',num2str(k-1),'[',num2str(0),'] = 0;\n']);
    fprintf(fidp0,'for(l%d=1;l%d<=%d;l%d++) {\n',k-1,k-1,lims_phi(k-1),k-1);
    if(deriv>=1)
        fprintf(fidp1,['phi_xt0_',num2str(k-1),'[',num2str(0),'] = 0;\n']);
        for j = 1:(nx*(1+nti)+1)
            fprintf(fidp1,['phi_xt1_',num2str(k-1),'[',num2str(j-1),'] = 0;\n']);
        end
        fprintf(fidp1,'for(l%d=1;l%d<=%d;l%d++) {\n',k-1,k-1,lims_phi(k-1),k-1);
    end
    if(deriv>=2)
        fprintf(fidp2,['phi_xt0_',num2str(k-1),'[',num2str(0),'] = 0;\n']);
        for j = 1:(nx*(1+nti)+1)
            for l = 1:(nx*(1+nti)+1)
                fprintf(fidp2,['phi_xt2_',num2str(k-1),'[',num2str(j-1+(l-1)*(nx*(1+nti)+1)),'] = 0;\n']);
            end
        end
        fprintf(fidp2,'for(l%d=1;l%d<=%d;l%d++) {\n',k-1,k-1,lims_phi(k-1),k-1);
    end
end
for k = size(phi_form,2):-1:1
    fprintf(fidp0,'%s\n',phichar{i,k}{~cellfun('isempty',phichar{i,k})});
    if(deriv>=1)
        fprintf(fidp1,'%s\n',dphidxtchar{i,k}{~cellfun('isempty',dphidxtchar{i,k})});
    end
    if(deriv>=2)
        fprintf(fidp2,'%s\n',d2phidxt2char{i,k}{~cellfun('isempty',d2phidxt2char{i,k})});
    end
    if(k>1)
        fprintf(fidp0,'}\n');
        if(deriv>=1)
            fprintf(fidp1,'%s\n',phichar{i,k}{~cellfun('isempty',phichar{i,k})});
            fprintf(fidp1,'}\n');
        end
        if(deriv>=2)
            fprintf(fidp2,'%s\n',phichar{i,k}{~cellfun('isempty',phichar{i,k})});
            fprintf(fidp2,'}\n');
        end
    end
end
fprintf(fidp0,'break;\n');
if(deriv>=1)
    fprintf(fidp1,'break;\n');
end
if(deriv>=2)
    fprintf(fidp2,'break;\n');
end
end
for i = 1:length(crhs)
    f{i} = subs([f_form;prhs{i}],clhs,crhs{i});
    fchar{i} = arrayfun(@ccode,formula(f{i}),'Uniform',0);
    fchar{i} = replace(fchar{i},'  t0 = ','');
    fchar{i} = replace(fchar{i},';','');
    for l = nf:-1:1
        fchar{i} = replace(fchar{i},['xl1_',num2str(l)],['x[',num2str(l-1),'+(l1-1)*',num2str(nf),']']);
    end
    for l = nx:-1:1
        fchar{i} = replace(fchar{i},['x',num2str(l)],['x[',num2str(l-1),']']);
    end
    for l = np:-1:1
        fchar{i} = replace(fchar{i},['p',num2str(l)],['p[',num2str(l-1),']']);
    end
    for l = 1:size(vars,2)
        fchar{i} = replace(fchar{i},[vars{1,l},'l1'],[vars{1,l},'[l1-1]']);
    end
    u{i} = crhs{i}(1:length(lb));
    uchar{i} = arrayfun(@ccode,formula(u{i}),'Uniform',0);
    uchar{i} = replace(uchar{i},'  t0 = ','');
    uchar{i} = replace(uchar{i},';','');
    for l = nx:-1:1
        uchar{i} = replace(uchar{i},['x',num2str(l)],['x[',num2str(l-1),']']);
    end
    for l = np:-1:1
        uchar{i} = replace(uchar{i},['p',num2str(l)],['p[',num2str(l-1),']']);
    end
    if(deriv>=1)
        dfdxp{i} = jacobian(f{i},[x(1:nx),reshape(xl(:,1:nf).',[],1).',p]);
        dfdxpchar{i} = arrayfun(@ccode,formula(dfdxp{i}),'Uniform',0);
        dfdxpchar{i} = replace(dfdxpchar{i},'  t0 = ','');
        dfdxpchar{i} = replace(dfdxpchar{i},';','');
        for l = nf:-1:1
            dfdxpchar{i} = replace(dfdxpchar{i},['xl1_',num2str(l)],['x[',num2str(l-1),'+(l1-1)*',num2str(nf),']']);
        end
        for l = nx:-1:1
            dfdxpchar{i} = replace(dfdxpchar{i},['x',num2str(l)],['x[',num2str(l-1),']']);
        end
        for l = np:-1:1
            dfdxpchar{i} = replace(dfdxpchar{i},['p',num2str(l)],['p[',num2str(l-1),']']);
        end
        for l = 1:size(vars,2)
            dfdxpchar{i} = replace(dfdxpchar{i},[vars{1,l},'l1'],[vars{1,l},'[l1-1]']);
        end
    else
        dfdxp{i} = [];
    end
    dudxp{i} = jacobian(u{i},x(1:nx));
    dudxpchar{i} = arrayfun(@ccode,formula(dudxp{i}),'Uniform',0);
    dudxpchar{i} = replace(dudxpchar{i},'  t0 = ','');
    dudxpchar{i} = replace(dudxpchar{i},';','');
    for l = nx:-1:1
        dudxpchar{i} = replace(dudxpchar{i},['x',num2str(l)],['x[',num2str(l-1),']']);
    end
    for l = np:-1:1
        dudxpchar{i} = replace(dudxpchar{i},['p',num2str(l)],['p[',num2str(l-1),']']);
    end
    if(deriv>=2)
        d2Hdxp2{i} = jacobian(dfdxp{i}.'*c(1:nx).',[x(1:nx),p]);
        for j = 1:nx
            d2fdxp2{i}{j} = diff(d2Hdxp2{i},c(j));
            d2fdxp2char{i}{j} = arrayfun(@ccode,formula(d2fdxp2{i}{j}),'Uniform',0);
            d2fdxp2char{i}{j} = replace(d2fdxp2char{i}{j},'  t0 = ','');
            d2fdxp2char{i}{j} = replace(d2fdxp2char{i}{j},';','');
            for l = nx:-1:1
                d2fdxp2char{i}{j} = replace(d2fdxp2char{i}{j},['x',num2str(l)],['x[',num2str(l-1),']']);
            end
            for l = np:-1:1
                d2fdxp2char{i}{j} = replace(d2fdxp2char{i}{j},['p',num2str(l)],['p[',num2str(l-1),']']);
            end
        end
        for j = 1:nx
            dndydtdunchar{i}{j} = ['f[',num2str(j-1),'] = ',fchar{i}{j},';'];
            for k = 1:(u_s-1)
                dndydtdunchar{i}{j+k*(nx+np)} = ['f[',num2str(j+k*(nx+np)-1),'] ='];
                for l = 1:nx
                    dndydtdunchar{i}{j+k*(nx+np)} = [dndydtdunchar{i}{j+k*(nx+np)},...
                        ' + (',dfdxpchar{i}{j,l},')*x[',num2str(l+k*(nx+np)-1),']'];
                end
                for l = 1:np
                    dndydtdunchar{i}{j+k*(nx+np)} = [dndydtdunchar{i}{j+k*(nx+np)},...
                        ' + (',dfdxpchar{i}{j,nx+nf*length(lims_phi)+l},')*x[',num2str(nx+l+k*(nx+np)-1),']'];
                end
                dndydtdunchar{i}{j+k*(nx+np)} = [dndydtdunchar{i}{j+k*(nx+np)},';'];
            end
        end
        for j = 1:nx
            for k = 0:(u_s-1)
                for l = 1:nx
                    dndadtdunchar{i}{j+k*(nx+np)+(l+k*(nx+np)-1)*(nx+np)*u_s} = ['f_x[',num2str(j+k*(nx+np)+(l+k*(nx+np)-1)*(nx+np)*u_s-1),'] ='];
                    dndadtdunchar{i}{j+k*(nx+np)+(l+k*(nx+np)-1)*(nx+np)*u_s} = [dndadtdunchar{i}{j+k*(nx+np)+(l+k*(nx+np)-1)*(nx+np)*u_s},...
                        ' + (',dfdxpchar{i}{j,l},')'];
                    dndadtdunchar{i}{j+k*(nx+np)+(l+k*(nx+np)-1)*(nx+np)*u_s} = [dndadtdunchar{i}{j+k*(nx+np)+(l+k*(nx+np)-1)*(nx+np)*u_s},';'];
                    if(k>0)
                        dndadtdunchar{i}{j+(l+k*(nx+np)-1)*(nx+np)*u_s} = ['f_x[',num2str(j+(l+k*(nx+np)-1)*(nx+np)*u_s-1),'] ='];
                        for m = 1:nx
                            dndadtdunchar{i}{j+(l+k*(nx+np)-1)*(nx+np)*u_s} = [dndadtdunchar{i}{j+(l+k*(nx+np)-1)*(nx+np)*u_s},...
                                ' + (',d2fdxp2char{i}{j}{l,m},')*x[',num2str(m+k*(nx+np)-1),']'];
                        end
                        for m = 1:np
                            dndadtdunchar{i}{j+(l+k*(nx+np)-1)*(nx+np)*u_s} = [dndadtdunchar{i}{j+(l+k*(nx+np)-1)*(nx+np)*u_s},...
                                ' + (',d2fdxp2char{i}{j}{l,nx+m},')*x[',num2str(nx+m+k*(nx+np)-1),']'];
                        end
                        dndadtdunchar{i}{j+(l+k*(nx+np)-1)*(nx+np)*u_s} = [dndadtdunchar{i}{j+(l+k*(nx+np)-1)*(nx+np)*u_s},';'];
                    end
                end
                for l = 1:np
                    dndadtdunchar{i}{j+k*(nx+np)+(nx+l+k*(nx+np)-1)*(nx+np)*u_s} = ['f_x[',num2str(j+k*(nx+np)+(nx+l+k*(nx+np)-1)*(nx+np)*u_s-1),'] ='];
                    dndadtdunchar{i}{j+k*(nx+np)+(nx+l+k*(nx+np)-1)*(nx+np)*u_s} = [dndadtdunchar{i}{j+k*(nx+np)+(nx+l+k*(nx+np)-1)*(nx+np)*u_s},...
                        ' + (',dfdxpchar{i}{j,nx+nf*length(lims_phi)+l},')'];
                    dndadtdunchar{i}{j+k*(nx+np)+(nx+l+k*(nx+np)-1)*(nx+np)*u_s} = [dndadtdunchar{i}{j+k*(nx+np)+(nx+l+k*(nx+np)-1)*(nx+np)*u_s},';'];
                    if(k>0)
                        dndadtdunchar{i}{j+(nx+l+k*(nx+np)-1)*(nx+np)*u_s} = ['f_x[',num2str(j+(nx+l+k*(nx+np)-1)*(nx+np)*u_s-1),'] ='];
                        for m = 1:nx
                            dndadtdunchar{i}{j+(nx+l+k*(nx+np)-1)*(nx+np)*u_s} = [dndadtdunchar{i}{j+(nx+l+k*(nx+np)-1)*(nx+np)*u_s},...
                                ' + (',d2fdxp2char{i}{j}{nx+l,m},')*x[',num2str(m+k*(nx+np)-1),']'];
                        end
                        for m = 1:np
                            dndadtdunchar{i}{j+(nx+l+k*(nx+np)-1)*(nx+np)*u_s} = [dndadtdunchar{i}{j+(nx+l+k*(nx+np)-1)*(nx+np)*u_s},...
                                ' + (',d2fdxp2char{i}{j}{nx+l,nx+m},')*x[',num2str(nx+m+k*(nx+np)-1),']'];
                        end
                        dndadtdunchar{i}{j+(nx+l+k*(nx+np)-1)*(nx+np)*u_s} = [dndadtdunchar{i}{j+(nx+l+k*(nx+np)-1)*(nx+np)*u_s},';'];
                    end
                end
            end
        end
    else
        d2Hdxp2{i} = [];
        for j = 1:nx
            d2fdxp2{i}{j} = [];
        end
    end
    for j = 1:(nf+ni*nu)
        if(~strcmp(fchar{i}{j},'0.0'))
            fchar{i}{j} = ['f[',num2str(j-1),'+(l1-1)*',num2str(nf),'] = ',fchar{i}{j},';'];
        else
            fchar{i}{j} = [];
        end
    end
    for j = 1:length(lb)
        if(~strcmp(uchar{i}{j},'0.0'))
            uchar{i}{j} = ['u[',num2str(j-1),'] = ',uchar{i}{j},';'];
        else
            uchar{i}{j} = [];
        end
    end
    if(deriv>=1)
        for j = 1:(nf+ni*nu)
            for l = 1:nx
                if(~strcmp(dfdxpchar{i}{j,l},'0.0'))
                    dfdxpchar{i}{j,l} = ['f_x[',num2str(j-1+(l-1)*(nx+np)),'+(l1-1)*',num2str(nf),'] += ',...
                        dfdxpchar{i}{j,l},';'];
                else
                    dfdxpchar{i}{j,l} = [];
                end
            end
            for l = 1:(nf*length(lims_phi))
                if(~strcmp(dfdxpchar{i}{j,nx+l},'0.0'))
                    dfdxpchar{i}{j,nx+l} = ['f_x[',num2str(j-1+(l-1)*(nx+np)),'+(l1-1)*',num2str(nf*(1+nx+np)),'] += ',...
                        dfdxpchar{i}{j,nx+l},';'];
                else
                    dfdxpchar{i}{j,nx+l} = [];
                end
            end
            for l = 1:np
                if(~strcmp(dfdxpchar{i}{j,nx+nf*length(lims_phi)+l},'0.0'))
                    dfdxpchar{i}{j,nx+nf*length(lims_phi)+l} = ['f_x[',num2str(j-1+(nx+l-1)*(nx+np)+(lims_f-1)*nf),'] = ',...
                        dfdxpchar{i}{j,nx+nf*length(lims_phi)+l},';'];
                else
                    dfdxpchar{i}{j,nx+nf*length(lims_phi)+l} = [];
                end
            end
        end
    end
    for j = 1:length(lb)
        for l = 1:nx
            if(~strcmp(dudxpchar{i}{j,l},'0.0'))
                dudxpchar{i}{j,l} = ['u_x[',num2str(j-1+(l-1)*length(lb)),'] = ',...
                    dudxpchar{i}{j,l},';'];
            else
                dudxpchar{i}{j,l} = [];
            end
        end
    end
    fprintf(fidy0,'case %d:\n',ci(i));
    fprintf(fidy0,'for(l1=1;l1<=%d;l1++) {\n',lims_f);
    fprintf(fidy0,'%s\n',fchar{i}{~cellfun('isempty',fchar{i})});
    fprintf(fidy0,'}\n');
    fprintf(fidy0,'if(u) {\n');
    fprintf(fidy0,'%s\n',uchar{i}{~cellfun('isempty',uchar{i})});
    fprintf(fidy0,'}\n');
    fprintf(fidy0,'break;\n');
    fprintf(fida0,'case %d:\n',ci(i));
    if(deriv>=1)
        fprintf(fida0,'for(l1=1;l1<=%d;l1++) {\n',lims_f);
        fprintf(fida0,'%s\n',dfdxpchar{i}{~cellfun('isempty',dfdxpchar{i})});
        fprintf(fida0,'}\n');
    end
    fprintf(fida0,'if(u_x) {\n');
    fprintf(fida0,'%s\n',dudxpchar{i}{~cellfun('isempty',dudxpchar{i})});
    fprintf(fida0,'}\n');
    fprintf(fida0,'break;\n');
    if(deriv>=2)
        fprintf(fidy1,'case %d:\n',ci(i));
        fprintf(fidy1,'for(l1=1;l1<=%d;l1++) {\n',lims_f);
        fprintf(fidy1,'%s\n',dndydtdunchar{i}{~cellfun('isempty',dndydtdunchar{i})});
        fprintf(fidy1,'}\n');
        fprintf(fidy1,'if(u) {\n');
        fprintf(fidy1,'%s\n',uchar{i}{~cellfun('isempty',uchar{i})});
        fprintf(fidy1,'}\n');
        fprintf(fidy1,'break;\n');
        fprintf(fida1,'case %d:\n',ci(i));
        fprintf(fida1,'for(l1=1;l1<=%d;l1++) {\n',lims_f);
        fprintf(fida1,'%s\n',dndadtdunchar{i}{~cellfun('isempty',dndadtdunchar{i})});
        fprintf(fida1,'}\n');
        fprintf(fida1,'if(u_x) {\n');
        fprintf(fida1,'%s\n',dudxpchar{i}{~cellfun('isempty',dudxpchar{i})});
        fprintf(fida1,'}\n');
        fprintf(fida1,'break;\n');
    end
end
h = h_form;
dhdx = jacobian(h,x(1:nx));
dhdu = jacobian(h,clhs);
hchar = arrayfun(@ccode,formula(h),'Uniform',0);
dhdxchar = arrayfun(@ccode,formula(dhdx),'Uniform',0);
dhduchar = arrayfun(@ccode,formula(dhdu),'Uniform',0);
hchar = replace(hchar,'  t0 = ','');
dhdxchar = replace(dhdxchar,'  t0 = ','');
dhduchar = replace(dhduchar,'  t0 = ','');
hchar = replace(hchar,';','');
dhdxchar = replace(dhdxchar,';','');
dhduchar = replace(dhduchar,';','');
for l = nx:-1:1
    hchar = replace(hchar,['x',num2str(l)],['x[',num2str(l-1),']']);
    dhdxchar = replace(dhdxchar,['x',num2str(l)],['x[',num2str(l-1),']']);
    dhduchar = replace(dhduchar,['x',num2str(l)],['x[',num2str(l-1),']']);
end
for l = nu:-1:1
    hchar = replace(hchar,['u',num2str(l)],['u[',num2str(l-1),']']);
    dhdxchar = replace(dhdxchar,['u',num2str(l)],['u[',num2str(l-1),']']);
    dhduchar = replace(dhduchar,['u',num2str(l)],['u[',num2str(l-1),']']);
end
for j = 1:nh
    if(~strcmp(hchar{j},'0.0'))
        hchar{j} = ['gh[',num2str(j-1),'] = ',hchar{j},';'];
    else
        hchar{j} = [];
    end
end
for j = 1:nh
    for l = 1:nx
        if(~strcmp(dhdxchar{j,l},'0.0')||~all(strcmp(dhduchar(j,:),'0.0')))
            dhdxchar{j,l} = ['gh_x[',num2str(j-1+(l-1)*nh),'] = ',dhdxchar{j,l}];
            for k = 1:nu
                dhdxchar{j,l} = [dhdxchar{j,l},'+(',dhduchar{j,k},')*u_x[',num2str(k-1+(l-1)*nu),']'];
            end
            dhdxchar{j,l} = [dhdxchar{j,l},';'];
        else
            dhdxchar{j,l} = [];
        end
    end
end
fprintf(fidgy,'%s\n',hchar{~cellfun('isempty',hchar)});
fprintf(fidga,'%s\n',dhdxchar{~cellfun('isempty',dhdxchar)});
fprintf(fidp0,'}\n}\n');
fprintf(fidy0,'}\n}\n');
fprintf(fida0,'}\n}\n');
if(deriv>=1)
    fprintf(fidp1,'}\n}\n');
end
if(deriv>=2)
    fprintf(fidp2,'}\n}\n');
    fprintf(fidy1,'}\n}\n');
    fprintf(fida1,'}\n}\n');
end
fprintf(fidgy,'}\n');
fprintf(fidga,'}\n');
fclose(fidp0);
fclose(fidy0);
fclose(fida0);
if(deriv>=1)
    fclose(fidp1);
end
if(deriv>=2)
    fclose(fidp2);
    fclose(fidy1);
    fclose(fida1);
end
fclose(fidgy);
fclose(fidga);
name = 'hessian';
eval(['mex -outdir ',name,' -output phixt0_',name,'_mex phixt0_mex.c ',name,'/phixt0_',name,'.c']);
eval(['mex -outdir ',name,' -output f_',name,'_mex f_mex.c ',name,'/f0_',name,'.c']);
eval(['mex -outdir ',name,' -output gh_',name,'_mex gh_mex.c ',name,'/f0_',name,'.c ',name,'/gh_',name,'.c']);
eval(['mex -outdir ',name,' -output ghx_',name,'_mex ghx_mex.c ',name,'/f0x_',name,'.c ',name,'/ghx_',name,'.c']);
eval(['mex -outdir ',name,' -output eventfun_',name,'_mex ../ode45_rp/eventfun_mex.c eventfun.c ',name,'/f0_',name,'.c ',name,'/gh_',name,'.c']);
eval(['mex -outdir ',name,' -output ode0_x_',name,'_mex ../ode45_rp/ode_mex.c ../ode45_rp/rp45.c ../ode45_rp/ntrp45.c odefun_x.c eventfun.c ',name,'/f0_',name,'.c ',name,'/gh_',name,'.c']);
if(deriv>=1)
    eval(['mex -outdir ',name,' -output phixt1_',name,'_mex phixt1_mex.c ',name,'/phixt1_',name,'.c']);
    eval(['mex -outdir ',name,' -output fx_',name,'_mex fx_mex.c ',name,'/f0x_',name,'.c']);
    eval(['mex -outdir ',name,' -output ode0_a_',name,'_mex ../ode45_rp/ode_mex.c ../ode45_rp/rp45.c ../ode45_rp/ntrp45.c odefun_a.c eventfun.c ',name,'/f0_',name,'.c ',name,'/gh_',name,'.c ',name,'/f0x_',name,'.c ../ode45_rp/binterp.c']);
end
if(deriv>=2)
    eval(['mex -outdir ',name,' -output phixt2_',name,'_mex phixt2_mex.c ',name,'/phixt2_',name,'.c']);
    eval(['mex -outdir ',name,' -output ode1_x_',name,'_mex ../ode45_rp/ode_mex.c ../ode45_rp/rp45.c ../ode45_rp/ntrp45.c odefun_x.c eventfun.c ',name,'/f1_',name,'.c ',name,'/gh_',name,'.c']);
    eval(['mex -outdir ',name,' -output ode1_a_',name,'_mex ../ode45_rp/ode_mex.c ../ode45_rp/rp45.c ../ode45_rp/ntrp45.c odefun_a.c eventfun.c ',name,'/f1_',name,'.c ',name,'/gh_',name,'.c ',name,'/f1x_',name,'.c ../ode45_rp/binterp.c']);
end
end