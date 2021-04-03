clc
close all
clear
profile on

% nh = 0;
% nx = 1;
% nh = 1;
nx = 2;
nh = 2;
if(nx==1)
    phi_form = str2sym({'xi1_1*x1^2';'t*x1*xi1_1'});
    f_form = str2sym({'u1'});
    x0 = str2sym({});
    lb = str2sym({});
    ub = str2sym({});
    gh_form = str2sym({});
    if(nh==0)
        h_form = str2sym({});
    else
        h_form = str2sym({'-x1+0.29'});
    end
    clhs = str2sym({'u1'});
    crhs = {str2sym({'-p1*x1'}),str2sym({'-p1*x1^2'}),str2sym({'-p1*x1^3'}),str2sym({'p1*x1^4'})};
    if(nh==0)
        nci = [];
    else
        nci = 4;
    end
    np = 1;
    ni = 0;
    nt = 3;
    nti = 2;
elseif(nx==2)
    phi_form = str2sym({'xi1_2*x1^2+xi1_1*x2^2';'t*x1*x2*xi2_1*xi2_2'});
    f_form = str2sym({'u1';'u2'});
    x0 = str2sym({});
    lb = str2sym({});
    ub = str2sym({});
    gh_form = str2sym({});
    if(nh==0)
        h_form = str2sym({});
    else
        h_form = str2sym({'x2-10*x1';'x2-20*x1'});
    end
    clhs = str2sym({'u1';'u2'});
    crhs = {str2sym({'-x1-2*p1*x1^2*x2';'x2-p1*x2^2*x1'}),...
        str2sym({'x1-3*p2*x1^2*x2';'2*x2-2*p2*x2^2*x1'}),...
        str2sym({'2*x1-p3*x1^2*x2';'-x2-3*p3*x2^2*x1'}),...
        str2sym({'x1+p1*x1^2*x2';'-x2+p1*x2^2*x1'})};
    if(nh==0)
        nci = [];
    else
        nci = [4,4];
    end
    np = 3;
    ni = 0;
    nt = 3;
    nti = 2;
end
phi_form = vpa(subs(phi_form));
f_form = vpa(subs(f_form));
x0 = double(subs(x0));
lb = double(subs(lb));
ub = double(subs(ub));
gh_form = vpa(subs(gh_form));
h_form = vpa(subs(h_form));
clhs = vpa(subs(clhs));
for i = 1:length(crhs)
    crhs{i} = vpa(subs(crhs{i}));
end
deriv = 2;

[nci,phi,dphidxt,d2phidxt2,f,dfdxp,d2Hdxp2,d2fdxp2,h,dhdx] = hessian_compile(nci,clhs,crhs,lb,ub,phi_form,f_form,gh_form,h_form,np,ni,nt,nti,deriv);
nphi = length(phi_form);

idx_x0 = length(x0)+(1:nx);
if(nx==1)
    av = [1,2,3];
    ti = [0.5,1.0];
	tsc = [0.8,1.4];
	tfc = 2.0;
    x0c = [x0;0.8];
    p0c = 1.0;
    sc = [1,1,1,1,1];
elseif(nx==2)
    av = [1,2,3];
    ti = [0.5,1.0];
	tsc = [0.8,1.4];
	tfc = 2.0;
    x0c = [x0;0.8;0.7];
    p0c = [1.0;1.0;1.0];
    sc = [1,1,1,1,1,1,1,1];
end
hessian_check(nci,1:nphi,av,idx_x0,ti,tsc,tfc,x0c,p0c,sc,deriv);

profile viewer