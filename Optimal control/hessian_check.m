function hessian_check(nci,phicase,av,idx_x0,ti,ts,tf,x0,p0,sc,deriv)
format long
lam = zeros(length(phicase),length(ts)+1+length(x0)+length(p0)+1);
pi0 = zeros(length(ts)+1+length(x0)+length(p0),1);
addpath('./hessian');
if(deriv<1)
    [J0,th0] = hessian_calc(nci,phicase,av,ti,0,ts,tf,x0,p0,lam,pi0);
    nphi = length(J0);
    nth = length(th0);
    rmpath('./hessian');
    for i = 1:nth
    th0i = cell2mat(th0(:,i));
    fprintf('Entry point:\n');
    disp(th0i);
    end
    for i = 1:nphi
    J0i = cell2mat(J0(:,i));
    fprintf('Cost:\n');
    disp(J0i);
    end
    return;
end
if(deriv>=2)
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,d2Jdtsdu,d2Jdtfdu,d2Jdx0du,d2Jdp0du,d2thdtsdu,d2thdtfdu,d2thdx0du,d2thdp0du,dtdulim,dxpdulim,dcqdulim,ddrdulim] =...
        hessian_calc(nci,phicase,av,ti,0,ts,tf,x0,p0,lam,pi0);
end
[J0,th0,dJdts0,dJdtf0,dJdx00,dJdp00,dthdts0,dthdtf0,dthdx00,dthdp00,tlim0,xplim0,cqlim0,drlim0] =...
    hessian_calc(nci,phicase,av,ti,0,ts,tf,x0,p0,lam,pi0);
nphi = length(J0);
nth = length(th0);
nt = length(ts)+1;
nx = length(x0);
np = length(p0);
nix = length(idx_x0);
draw_plot = true;
if(draw_plot==true)
    hessian_plot(nci,phicase,av,ti,0,ts,tf,x0,p0,lam,pi0,[]);
end
pd = sc*1e-6;
tlimpts = cell(nt-1,1);
tlimptf = cell(1,1);
tlimpx0 = cell(nix,1);
tlimpp0 = cell(np,1);
xplimpts = cell(nt-1,1);
xplimptf = cell(1,1);
xplimpx0 = cell(nix,1);
xplimpp0 = cell(np,1);
thpts = cell(nt-1,nth);
thptf = cell(1,nth);
thpx0 = cell(nix,nth);
thpp0 = cell(np,nth);
dthdtspts = cell(nt-1,nth);
dthdtsptf = cell(1,nth);
dthdtspx0 = cell(nix,nth);
dthdtspp0 = cell(np,nth);
dthdtfpts = cell(nt-1,nth);
dthdtfptf = cell(1,nth);
dthdtfpx0 = cell(nix,nth);
dthdtfpp0 = cell(np,nth);
dthdx0pts = cell(nt-1,nth);
dthdx0ptf = cell(1,nth);
dthdx0px0 = cell(nix,nth);
dthdx0pp0 = cell(np,nth);
dthdp0pts = cell(nt-1,nth);
dthdp0ptf = cell(1,nth);
dthdp0px0 = cell(nix,nth);
dthdp0pp0 = cell(np,nth);
drlimpts = cell(nt-1,nth);
drlimptf = cell(1,nth);
drlimpx0 = cell(nix,nth);
drlimpp0 = cell(np,nth);
Jpts = cell(nt-1,nphi);
Jptf = cell(1,nphi);
Jpx0 = cell(nix,nphi);
Jpp0 = cell(np,nphi);
dJdtspts = cell(nt-1,nphi);
dJdtsptf = cell(1,nphi);
dJdtspx0 = cell(nix,nphi);
dJdtspp0 = cell(np,nphi);
dJdtfpts = cell(nt-1,nphi);
dJdtfptf = cell(1,nphi);
dJdtfpx0 = cell(nix,nphi);
dJdtfpp0 = cell(np,nphi);
dJdx0pts = cell(nt-1,nphi);
dJdx0ptf = cell(1,nphi);
dJdx0px0 = cell(nix,nphi);
dJdx0pp0 = cell(np,nphi);
dJdp0pts = cell(nt-1,nphi);
dJdp0ptf = cell(1,nphi);
dJdp0px0 = cell(nix,nphi);
dJdp0pp0 = cell(np,nphi);
cqlimpts = cell(nt-1,nphi);
cqlimptf = cell(1,nphi);
cqlimpx0 = cell(nix,nphi);
cqlimpp0 = cell(np,nphi);
for k = 1:(nt-1)
    tsp = ts;
    tsp(k) = tsp(k)+pd(k);
    [Jpts(k,:),thpts(k,:),dJdtspts(k,:),dJdtfpts(k,:),dJdx0pts(k,:),dJdp0pts(k,:),dthdtspts(k,:),dthdtfpts(k,:),dthdx0pts(k,:),dthdp0pts(k,:),tlimpts(k,:),xplimpts(k,:),cqlimpts(k,:),drlimpts(k,:,:)] =...
        hessian_calc(nci,phicase,av,ti,0,tsp,tf,x0,p0,lam,pi0);
end
tfp = tf+pd(nt);
[Jptf(:,:),thptf(:,:),dJdtsptf(:,:),dJdtfptf(:,:),dJdx0ptf(:,:),dJdp0ptf(:,:),dthdtsptf(:,:),dthdtfptf(:,:),dthdx0ptf(:,:),dthdp0ptf(:,:),tlimptf(:,:),xplimptf(:,:),cqlimptf(:,:),drlimptf(:,:,:)] =...
    hessian_calc(nci,phicase,av,ti,0,ts,tfp,x0,p0,lam,pi0);
for k = 1:nix
    x0p = x0;
    x0p(idx_x0(k)) = x0p(idx_x0(k))+pd(nt+k);
    [Jpx0(k,:),thpx0(k,:),dJdtspx0(k,:),dJdtfpx0(k,:),dJdx0px0(k,:),dJdp0px0(k,:),dthdtspx0(k,:),dthdtfpx0(k,:),dthdx0px0(k,:),dthdp0px0(k,:),tlimpx0(k,:),xplimpx0(k,:),cqlimpx0(k,:),drlimpx0(k,:,:)] =...
        hessian_calc(nci,phicase,av,ti,0,ts,tf,x0p,p0,lam,pi0);
end
for k = 1:np
    p0p = p0;
    p0p(k) = p0p(k)+pd(nt+nix+k);
    [Jpp0(k,:),thpp0(k,:),dJdtspp0(k,:),dJdtfpp0(k,:),dJdx0pp0(k,:),dJdp0pp0(k,:),dthdtspp0(k,:),dthdtfpp0(k,:),dthdx0pp0(k,:),dthdp0pp0(k,:),tlimpp0(k,:),xplimpp0(k,:),cqlimpp0(k,:),drlimpp0(k,:,:)] =...
        hessian_calc(nci,phicase,av,ti,0,ts,tf,x0,p0p,lam,pi0);
end
rmpath('./hessian');
tlim0 = cell2mat(cellfun(@(e)shiftdim(e,-1),tlim0,'UniformOutput',false));
tlimpts = cell2mat(cellfun(@(e)shiftdim(e,-1),tlimpts,'UniformOutput',false));
tlimptf = cell2mat(cellfun(@(e)shiftdim(e,-1),tlimptf,'UniformOutput',false));
tlimpx0 = cell2mat(cellfun(@(e)shiftdim(e,-1),tlimpx0,'UniformOutput',false));
tlimpp0 = cell2mat(cellfun(@(e)shiftdim(e,-1),tlimpp0,'UniformOutput',false));
xplim0 = cell2mat(cellfun(@(e)shiftdim(e,-1),xplim0,'UniformOutput',false));
xplimpts = cell2mat(cellfun(@(e)shiftdim(e,-1),xplimpts,'UniformOutput',false));
xplimptf = cell2mat(cellfun(@(e)shiftdim(e,-1),xplimptf,'UniformOutput',false));
xplimpx0 = cell2mat(cellfun(@(e)shiftdim(e,-1),xplimpx0,'UniformOutput',false));
xplimpp0 = cell2mat(cellfun(@(e)shiftdim(e,-1),xplimpp0,'UniformOutput',false));
if(deriv>=2)
    dtdulim = cell2mat(dtdulim);
    dxpdulim = cell2mat(dxpdulim);
    fprintf('Derivatives:\n');
    fprintf('Switching times:\n');
    dtdulim = dtdulim(:,[1:nt,nt+idx_x0,nt+nx+(1:np)],:);
    dtlimpts = permute(reshape(tlimpts,[],size(tlim0,2),size(tlim0,3))-tlim0,[2,1,3]);
    dtlimptf = permute(reshape(tlimptf,[],size(tlim0,2),size(tlim0,3))-tlim0,[2,1,3]);
    dtlimpx0 = permute(reshape(tlimpx0,[],size(tlim0,2),size(tlim0,3))-tlim0,[2,1,3]);
    dtlimpp0 = permute(reshape(tlimpp0,[],size(tlim0,2),size(tlim0,3))-tlim0,[2,1,3]);
    dtdulimn = cat(2,dtlimpts,dtlimptf,dtlimpx0,dtlimpp0)./pd;
    disp((dtdulimn-dtdulim).*sc);
    fprintf('States:\n');
    dxpdulim = dxpdulim(:,[1:nt,nt+idx_x0,nt+nx+(1:np)],:);
    dxplimpts = permute(reshape(xplimpts,[],size(xplim0,2),size(xplim0,3))-xplim0,[2,1,3]);
    dxplimptf = permute(reshape(xplimptf,[],size(xplim0,2),size(xplim0,3))-xplim0,[2,1,3]);
    dxplimpx0 = permute(reshape(xplimpx0,[],size(xplim0,2),size(xplim0,3))-xplim0,[2,1,3]);
    dxplimpp0 = permute(reshape(xplimpp0,[],size(xplim0,2),size(xplim0,3))-xplim0,[2,1,3]);
    dxpdulimn = cat(2,dxplimpts,dxplimptf,dxplimpx0,dxplimpp0)./pd;
    disp((dxpdulimn-dxpdulim).*sc);
end
for i = 1:nth
th0i = cell2mat(th0(:,i));
thptsi = cell2mat(thpts(:,i));
thptfi = cell2mat(thptf(:,i));
thpx0i = cell2mat(thpx0(:,i));
thpp0i = cell2mat(thpp0(:,i));
dthdts0i = cell2mat(dthdts0(:,i));
dthdtsptsi = cell2mat(dthdtspts(:,i));
dthdtsptfi = cell2mat(dthdtsptf(:,i));
dthdtspx0i = cell2mat(dthdtspx0(:,i));
dthdtspp0i = cell2mat(dthdtspp0(:,i));
dthdtf0i = cell2mat(dthdtf0(:,i));
dthdtfptsi = cell2mat(dthdtfpts(:,i));
dthdtfptfi = cell2mat(dthdtfptf(:,i));
dthdtfpx0i = cell2mat(dthdtfpx0(:,i));
dthdtfpp0i = cell2mat(dthdtfpp0(:,i));
dthdx00i = cell2mat(dthdx00(:,i));
dthdx0ptsi = cell2mat(dthdx0pts(:,i));
dthdx0ptfi = cell2mat(dthdx0ptf(:,i));
dthdx0px0i = cell2mat(dthdx0px0(:,i));
dthdx0pp0i = cell2mat(dthdx0pp0(:,i));
dthdp00i = cell2mat(dthdp00(:,i));
dthdp0ptsi = cell2mat(dthdp0pts(:,i));
dthdp0ptfi = cell2mat(dthdp0ptf(:,i));
dthdp0px0i = cell2mat(dthdp0px0(:,i));
dthdp0pp0i = cell2mat(dthdp0pp0(:,i));
if(deriv>=2)
    d2thdtsdui = cell2mat(d2thdtsdu(:,i));
    d2thdtfdui = cell2mat(d2thdtfdu(:,i));
    d2thdx0dui = cell2mat(d2thdx0du(:,i));
    d2thdp0dui = cell2mat(d2thdp0du(:,i));
end
drlim0i = cell2mat(cellfun(@(e)shiftdim(e,-1),drlim0(:,i),'UniformOutput',false));
drlimptsi = cell2mat(cellfun(@(e)shiftdim(e,-1),drlimpts(:,i),'UniformOutput',false));
drlimptfi = cell2mat(cellfun(@(e)shiftdim(e,-1),drlimptf(:,i),'UniformOutput',false));
drlimpx0i = cell2mat(cellfun(@(e)shiftdim(e,-1),drlimpx0(:,i),'UniformOutput',false));
drlimpp0i = cell2mat(cellfun(@(e)shiftdim(e,-1),drlimpp0(:,i),'UniformOutput',false));
fprintf('Entry point:\n');
disp(th0i);
fprintf('Gradient of entry point:\n');
dthdu = [dthdts0i,dthdtf0i,dthdx00i(idx_x0),dthdp00i];
disp(dthdu)
dthdun = ([thptsi;thptfi;thpx0i;thpp0i]-th0i)'./pd;
disp((dthdun-dthdu)./dthdu);
if(deriv>=2)
    fprintf('Hessian of entry point:\n');
    d2thdu2 = [d2thdtsdui(:,[1:nt,nt+idx_x0,nt+nx+(1:np)]);...
        d2thdtfdui(:,[1:nt,nt+idx_x0,nt+nx+(1:np)]);...
        d2thdx0dui(idx_x0,[1:nt,nt+idx_x0,nt+nx+(1:np)]);...
        d2thdp0dui(:,[1:nt,nt+idx_x0,nt+nx+(1:np)])];
    disp(d2thdu2)
    d2thdu2n = [[dthdtsptsi;dthdtsptfi;dthdtspx0i;dthdtspp0i]-dthdts0i,...
        [dthdtfptsi;dthdtfptfi;dthdtfpx0i;dthdtfpp0i]-dthdtf0i,...
        [dthdx0ptsi(:,idx_x0);dthdx0ptfi(:,idx_x0);dthdx0px0i(:,idx_x0);dthdx0pp0i(:,idx_x0)]-dthdx00i(idx_x0),...
        [dthdp0ptsi;dthdp0ptfi;dthdp0px0i;dthdp0pp0i]-dthdp00i]'...
        ./pd;
    disp((d2thdu2n-d2thdu2)./d2thdu2);
    ddrdulimi = cell2mat(ddrdulim(:,i));
    fprintf('Derivatives:\n');
    fprintf('Adjoint variables for entry point:\n');
    ddrdulimi = ddrdulimi(:,[1:nt,nt+idx_x0,nt+nx+(1:np)],:);
    ddrlimpts = permute(reshape(drlimptsi,[],size(drlim0i,2),size(drlim0i,3))-drlim0i,[2,1,3]);
    ddrlimptf = permute(reshape(drlimptfi,[],size(drlim0i,2),size(drlim0i,3))-drlim0i,[2,1,3]);
    ddrlimpx0 = permute(reshape(drlimpx0i,[],size(drlim0i,2),size(drlim0i,3))-drlim0i,[2,1,3]);
    ddrlimpp0 = permute(reshape(drlimpp0i,[],size(drlim0i,2),size(drlim0i,3))-drlim0i,[2,1,3]);
    ddrdulimn = cat(2,ddrlimpts,ddrlimptf,ddrlimpx0,ddrlimpp0)./pd;
    disp((ddrdulimn-ddrdulimi).*sc/tf);
end
end
for i = 1:nphi
J0i = cell2mat(J0(:,i));
Jptsi = cell2mat(Jpts(:,i));
Jptfi = cell2mat(Jptf(:,i));
Jpx0i = cell2mat(Jpx0(:,i));
Jpp0i = cell2mat(Jpp0(:,i));
dJdts0i = cell2mat(dJdts0(:,i));
dJdtsptsi = cell2mat(dJdtspts(:,i));
dJdtsptfi = cell2mat(dJdtsptf(:,i));
dJdtspx0i = cell2mat(dJdtspx0(:,i));
dJdtspp0i = cell2mat(dJdtspp0(:,i));
dJdtf0i = cell2mat(dJdtf0(:,i));
dJdtfptsi = cell2mat(dJdtfpts(:,i));
dJdtfptfi = cell2mat(dJdtfptf(:,i));
dJdtfpx0i = cell2mat(dJdtfpx0(:,i));
dJdtfpp0i = cell2mat(dJdtfpp0(:,i));
dJdx00i = cell2mat(dJdx00(:,i));
dJdx0ptsi = cell2mat(dJdx0pts(:,i));
dJdx0ptfi = cell2mat(dJdx0ptf(:,i));
dJdx0px0i = cell2mat(dJdx0px0(:,i));
dJdx0pp0i = cell2mat(dJdx0pp0(:,i));
dJdp00i = cell2mat(dJdp00(:,i));
dJdp0ptsi = cell2mat(dJdp0pts(:,i));
dJdp0ptfi = cell2mat(dJdp0ptf(:,i));
dJdp0px0i = cell2mat(dJdp0px0(:,i));
dJdp0pp0i = cell2mat(dJdp0pp0(:,i));
if(deriv>=2)
    d2Jdtsdui = cell2mat(d2Jdtsdu(:,i));
    d2Jdtfdui = cell2mat(d2Jdtfdu(:,i));
    d2Jdx0dui = cell2mat(d2Jdx0du(:,i));
    d2Jdp0dui = cell2mat(d2Jdp0du(:,i));
end
cqlim0i = cell2mat(cellfun(@(e)shiftdim(e,-1),cqlim0(:,i),'UniformOutput',false));
cqlimptsi = cell2mat(cellfun(@(e)shiftdim(e,-1),cqlimpts(:,i),'UniformOutput',false));
cqlimptfi = cell2mat(cellfun(@(e)shiftdim(e,-1),cqlimptf(:,i),'UniformOutput',false));
cqlimpx0i = cell2mat(cellfun(@(e)shiftdim(e,-1),cqlimpx0(:,i),'UniformOutput',false));
cqlimpp0i = cell2mat(cellfun(@(e)shiftdim(e,-1),cqlimpp0(:,i),'UniformOutput',false));
fprintf('Cost:\n');
disp(J0i);
fprintf('Gradient of cost:\n');
dJdu = [dJdts0i,dJdtf0i,dJdx00i(idx_x0),dJdp00i];
disp(dJdu)
dJdun = ([Jptsi;Jptfi;Jpx0i;Jpp0i]-J0i)'./pd;
disp((dJdun-dJdu)./dJdu);
if(deriv>=2)
    fprintf('Hessian of cost:\n');
    d2Jdu2 = [d2Jdtsdui(:,[1:nt,nt+idx_x0,nt+nx+(1:np)]);...
        d2Jdtfdui(:,[1:nt,nt+idx_x0,nt+nx+(1:np)]);...
        d2Jdx0dui(idx_x0,[1:nt,nt+idx_x0,nt+nx+(1:np)]);...
        d2Jdp0dui(:,[1:nt,nt+idx_x0,nt+nx+(1:np)])];
    disp(d2Jdu2)
    d2Jdu2n = [[dJdtsptsi;dJdtsptfi;dJdtspx0i;dJdtspp0i]-dJdts0i,...
        [dJdtfptsi;dJdtfptfi;dJdtfpx0i;dJdtfpp0i]-dJdtf0i,...
        [dJdx0ptsi(:,idx_x0);dJdx0ptfi(:,idx_x0);dJdx0px0i(:,idx_x0);dJdx0pp0i(:,idx_x0)]-dJdx00i(idx_x0),...
        [dJdp0ptsi;dJdp0ptfi;dJdp0px0i;dJdp0pp0i]-dJdp00i]'...
        ./pd;
    disp((d2Jdu2n-d2Jdu2)./d2Jdu2);
    dcqdulimi = cell2mat(dcqdulim(:,i));
    fprintf('Derivatives:\n');
    fprintf('Adjoint variables for cost:\n');
    dcqdulimi = dcqdulimi(:,[1:nt,nt+idx_x0,nt+nx+(1:np)],:);
    dcqlimpts = permute(reshape(cqlimptsi,[],size(cqlim0i,2),size(cqlim0i,3))-cqlim0i,[2,1,3]);
    dcqlimptf = permute(reshape(cqlimptfi,[],size(cqlim0i,2),size(cqlim0i,3))-cqlim0i,[2,1,3]);
    dcqlimpx0 = permute(reshape(cqlimpx0i,[],size(cqlim0i,2),size(cqlim0i,3))-cqlim0i,[2,1,3]);
    dcqlimpp0 = permute(reshape(cqlimpp0i,[],size(cqlim0i,2),size(cqlim0i,3))-cqlim0i,[2,1,3]);
    dcqdulimn = cat(2,dcqlimpts,dcqlimptf,dcqlimpx0,dcqlimpp0)./pd;
    disp((dcqdulimn-dcqdulimi).*sc);
end
end
end