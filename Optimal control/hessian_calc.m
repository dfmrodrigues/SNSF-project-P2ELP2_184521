function [J,th,dJdts,dJdtf,dJdx0,dJdp0,dthdts,dthdtf,dthdx0,dthdp0,tlim,xplim,cqlim,drlim,tout,xpout,cqout,drout,uout,...
    d2Jdtsdu,d2Jdtfdu,d2Jdx0du,d2Jdp0du,d2thdtsdu,d2thdtfdu,d2thdx0du,d2thdp0du,dtdulim,dxpdulim,dcqdulim,ddrdulim] = hessian_calc(nci,phicase,av,ti,t0,ts,tf,x0,p0,lam,pi0)
tv = [ts,tf];
nt = length(tv);
[tiv,iiv] = sort([tv,ti]);
nti = length(tiv);
nx = length(x0);
nu = size(av,1);
np = length(p0);
nh = length(nci);
if(nargout>19)
    second = 1;
else
    second = 0;
end
if(nargout>10||nargout<=2)
    fast = 0;
else
    fast = 1;
end
if(nargout>2)
    first = 1;
else
    first = 0;
end
iev = [];
civ = [];
tlim = t0;
xlim = [];
if(second)
    dtdulim = [];
    dxdulim = [];
end
if(second)
    dpdu = [zeros(np,nt+nx),eye(np)];
    y0 = reshape([x0,[zeros(nx,nt),eye(nx),zeros(nx,np)];p0,dpdu],[],1);
else
    y0 = [x0;p0];
end
if(fast)
    ip = 5000;
else
    ip = 50000;
end
deltat = tf-t0;
tout = zeros(1,nti*ip+1);
yout = zeros(length(y0),length(tout));
uout = zeros(nu,length(tout));
idx_tout = 1;
tout(idx_tout) = t0;
yout(:,idx_tout) = y0;
for k = 1:nti
    tk = tout((k-1)*ip+1);
    if(tk-tiv(k)>-1e-12)
        tiv(k) = tk+1e-12;
        if(iiv(k)<=nt)
            tv(iiv(k)) = tiv(k);
        end
    end
    tout(((k-1)*ip+1):(k*ip+1)) = tk:(tiv(k)-tk)/ip:tiv(k);
end
if(second)
    dtdu = [tv==t0,zeros(1,nx+np)];
end
for k = 1:nt
    cases = av(k);
    civ = cat(2,civ,cases);
    while(true)
        [gh0,f0,~] = gh_hessian_mex([x0;p0],nu,nh,cases);
        ghx0 = ghx_hessian_mex([x0;p0],nu,nh,cases);
        ie = find(or(gh0>1e-8,and(gh0>-1e-8,ghx0*f0>1e-8)));
        if(isempty(ie))
            break;
        end
        cases = nci(ie);
        civ(:,end) = cases;
    end
    while(t0<tv(k))
        tspan = tout(idx_tout(end):(find(iiv==k,1,'first')*ip+1));
        if(second)
            dxpdu = reshape(y0((nx+np+1):end),nx+np,[]);
            dxpdu = dxpdu-f_hessian_mex([x0;p0],nu,cases)*dtdu;
            y0((nx+np+1):end) = reshape(dxpdu,[],1);
            if(length(tspan)>1&&all(diff(tspan)>0))
                [toutk,youtk,~,~,ie] = ode45_rp(@ode1_x_hessian_mex,tspan,y0,1e-9,1e-12,{@eventfun_hessian_mex,nh},nu,cases);
            else
                toutk = tv(k);
                youtk = y0;
                ie = 0;
            end
        else
            if(length(tspan)>1&&all(diff(tspan)>0))
                [toutk,youtk,~,~,ie] = ode45_rp(@ode0_x_hessian_mex,tspan,y0,1e-9,1e-12,{@eventfun_hessian_mex,nh},nu,cases);
            else
                toutk = tv(k);
                youtk = y0;
                ie = 0;
            end
        end
        if(~fast)
            uoutk = zeros(nu,length(toutk));
            for j = 1:length(toutk)
                [~,uoutk(:,j)] = f_hessian_mex(youtk(1:(nx+np),j),nu,cases);
            end
        end
        t0 = tout(idx_tout(end));
        y0 = yout(:,idx_tout(end));
        idx_tout = idx_tout(end)+(0:(length(toutk)-1));
        tout(idx_tout) = toutk;
        yout(:,idx_tout) = youtk;
        if(~fast)
            uout(:,idx_tout) = uoutk;
        end
        tout(idx_tout(1)) = t0;
        yout(:,idx_tout(1)) = y0;
        t0 = toutk(end);
        y0 = youtk(:,end);
        if(t0==tv(k))
            ie = 0;
        end
        x0 = y0(1:nx);
        if(ie)
            ie = ie(1);
            dhadxp = (ie==(1:nh))*ghx_hessian_mex([x0;p0],nu,nh,cases);
            dhadt = dhadxp*f_hessian_mex([x0;p0],nu,cases);
        end
        if(second)
            dxpdu = reshape(y0((nx+np+1):end),nx+np,[]);
            if(ie&&prod(civ(:,end)==nci(ie)))
                dtdu = zeros(1,nt+nx+np);
            elseif(ie)
                dtdu = -dhadxp*dxpdu/dhadt;
            else
                dtdu = [tv==t0,zeros(1,nx+np)];
            end
            dxpdu = dxpdu+f_hessian_mex([x0;p0],nu,cases)*dtdu;
            y0((nx+np+1):end) = reshape(dxpdu,[],1);
        end
        if(ie&&prod(civ(:,end)==nci(ie)))
            continue;
        end
        iev = cat(2,iev,ie);
        if(ie)
            cases = nci(ie);
            civ = cat(2,civ,cases);
        end
        tlim = cat(2,tlim,t0);
        xlim = cat(2,xlim,x0);
        if(second)
            dtdulim = cat(3,dtdulim,dtdu);
            dxdulim = cat(3,dxdulim,dxpdu(1:nx,:));
        end
    end
end
[tlim,ilim] = sort([tlim,ti]);
xpout = yout(1:(nx+np),:);
xpi = reshape(cell2mat(arrayfun(@(t)binterp_mex(tout,xpout,t),ti,'UniformOutput',false)),nx+np,[]);
xi = reshape(xpi(1:nx,:),[],1);
xlim = cat(2,xlim,xpi(1:nx,:));
xlim = xlim(:,ilim(2:end)-1);
if(second)
    dtdulim = cat(3,dtdulim,zeros(1,nt+nx+np,length(ti)));
    dtdulim = dtdulim(:,:,ilim(2:end)-1);
    dxpduout = yout((nx+np+1):end,:);
    dxpdui = reshape(cell2mat(arrayfun(@(t)binterp_mex(tout,dxpduout,t),ti,'UniformOutput',false)),nx+np,nt+nx+np,[]);
    dxdulim = cat(3,dxdulim,dxpdui(1:nx,:,:));
    dxdulim = dxdulim(:,:,ilim(2:end)-1);
end
iev = cat(2,iev,zeros(1,length(ti)));
iev = iev(:,ilim(2:end)-1);
civ = cat(2,civ,zeros(1,length(ti)));
civ = civ(:,ilim(2:end)-1);
etv = zeros(1,length(tlim)-1);
for j = 1:(length(tlim)-1)
    etv(j) = j-1+find(civ(j:end),1,'first');
end
civ = civ(etv);
ith = find(iev);
nphi = length(phicase);
nth = length(ith);
J = cell(1,nphi);
th = cell(1,nth);
tf = tlim(end);
xf = xlim(:,end);
for i = 1:nphi
    J{i} = phixt0_hessian_mex(tf,xf,xi,phicase(i))+lam(i,1)+lam(i,2:end)*([ts';tf;x0;p0]-pi0);
end
for i = 1:nth
    th{i} = tlim(ith(i)+1);
end
if(~first)
    return;
end
dJdxt = cell(1,nphi);
dJdth = cell(1,nphi);
dJdts = cell(1,nphi);
dJdtf = cell(1,nphi);
dJdx0 = cell(1,nphi);
dJdp0 = cell(1,nphi);
cqout = cell(1,nphi);
cqlim = cell(1,nphi);
clim = cell(1,nphi);
qlim = cell(1,nphi);
af = cell(1,nphi);
qtf = cell(1,nphi);
aout = cell(1,nphi);
dthdth = cell(1,nth);
dthdts = cell(1,nth);
dthdtf = cell(1,nth);
dthdx0 = cell(1,nth);
dthdp0 = cell(1,nth);
drout = cell(1,nth);
drlim = cell(1,nth);
dlim = cell(1,nth);
rlim = cell(1,nth);
bf = cell(1,nth);
rtf = cell(1,nth);
bout = cell(1,nth);
if(second)
    d2Jdxt2 = cell(1,nphi);
    d2Jdthdu = cell(1,nphi);
    d2Jdtsdu = cell(1,nphi);
    d2Jdtfdu = cell(1,nphi);
    d2Jdx0du = cell(1,nphi);
    d2Jdp0du = cell(1,nphi);
    dcqdulim = cell(1,nphi);
    dcdulim = cell(1,nphi);
    dqdulim = cell(1,nphi);
    d2thdthdu = cell(1,nth);
    d2thdtsdu = cell(1,nth);
    d2thdtfdu = cell(1,nth);
    d2thdx0du = cell(1,nth);
    d2thdp0du = cell(1,nth);
    ddrdulim = cell(1,nth);
    dddulim = cell(1,nth);
    drdulim = cell(1,nth);
end
if(second)
    dtdu = [tv==tf,zeros(1,nx+np)];
end
cases = civ(end);
for i = 1:nphi
    clim{i} = [];
    qlim{i} = [];
    if(second)
        dcdulim{i} = [];
        dqdulim{i} = [];
    end
    dJdxt{i} = phixt1_hessian_mex(tf,xf,xi,phicase(i));
    if(second)
        d2Jdxt2{i} = phixt2_hessian_mex(tf,xf,xi,phicase(i));
        dxidu = reshape(permute(dxdulim(:,:,ismember(tlim(2:end),ti)),[2,1,3]),nt+nx+np,[])';
        af{i} = reshape([dJdxt{i}(2:(nx+1))',d2Jdxt2{i}(2:(nx+1),:)*[[zeros(1,nt-1),eye(1,1+nx+np)];dxdulim(:,:,end);dxidu];zeros(np,1),zeros(np,nt+nx+np)],[],1);
    else
        af{i} = [dJdxt{i}(2:(nx+1))';zeros(np,1)];
        aout{i} = zeros(length(af{i}),length(tout));
        aout{i}(:,end) = af{i};
    end
    qtf{i} = dJdxt{i}(1);
end
for i = 1:nth
    dlim{i} = [];
    rlim{i} = [];
    if(second)
        dddulim{i} = [];
        drdulim{i} = [];
    end
    if(second)
        bf{i} = reshape([zeros(nx+np,1),zeros(nx+np,nt+nx+np)],[],1);
    else
        bf{i} = zeros(nx+np,1);
        bout{i} = zeros(length(bf{i}),length(tout));
        bout{i}(:,end) = bf{i};
    end
    rtf{i} = 0;
end
for k = (length(tlim)-1):-1:1
    t0 = tlim(k);
    idx_tout = and(tout>=t0,tout<=tf);
    tspan = fliplr(tout(idx_tout));
    if(k==length(tlim)-1)
        for i = 1:nphi
            clim{i} = cat(2,af{i}(1:nx),clim{i});
            qlim{i} = cat(2,af{i}((nx+1):(nx+np)),qlim{i});
            if(second)
                dcqdu{i} = reshape(af{i}((nx+np+1):end),nx+np,[]);
                dcdulim{i} = cat(3,dcqdu{i}(1:nx,:),dcdulim{i});
                dqdulim{i} = cat(3,dcqdu{i}((nx+1):(nx+np),:),dqdulim{i});
            end
        end
        for i = 1:nth
            dlim{i} = cat(2,bf{i}(1:nx),dlim{i});
            rlim{i} = cat(2,bf{i}((nx+1):(nx+np)),rlim{i});
            if(second)
                ddrdu{i} = reshape(bf{i}((nx+np+1):end),nx+np,[]);
                dddulim{i} = cat(3,ddrdu{i}(1:nx,:),dddulim{i});
                drdulim{i} = cat(3,ddrdu{i}((nx+1):(nx+np),:),drdulim{i});
            end
        end
    end
    for i = 1:nphi
        af{i}(1:nx) = af{i}(1:nx)+kron(ti==tf,eye(nx))*dJdxt{i}((nx+2):end)';
        if(second)
            dcqdu{i} = dcqdu{i}+fx_hessian_mex([xf;p0],nu,cases)'*af{i}(1:(nx+np))*dtdu+...
                [kron(ti==tf,eye(nx))*d2Jdxt2{i}((nx+2):end,:)*[[zeros(1,nt-1),eye(1,1+nx+np)];dxdulim(:,:,end);dxidu];zeros(np,nt+nx+np)];
            af{i}((nx+np+1):end) = reshape(dcqdu{i},[],1);
            if(length(tspan)>1&&all(-diff(tspan)>0))
                [~,aoutk{i}] = ode45_rp(@ode1_a_hessian_mex,tspan,[af{i};qtf{i}],1e-9,1e-12,[],tout,yout,deltat,nu,cases);
            else
                aoutk{i} = [af{i};qtf{i}];
            end
        else
            if(length(tspan)>1&&all(-diff(tspan)>0))
                [~,aoutk{i}] = ode45_rp(@ode0_a_hessian_mex,tspan,[af{i};qtf{i}],1e-9,1e-12,[],tout,yout,deltat,nu,cases);
            else
                aoutk{i} = [af{i};qtf{i}];
            end
        end
        af{i} = aoutk{i}(1:end-1,end);
        qtf{i} = aoutk{i}(end,end);
        if(~fast)
            aoutk{i} = aoutk{i}(1:end-1,:);
            aout{i}(:,idx_tout) = fliplr(aoutk{i});
        end
        if(second)
            dcqdu{i} = reshape(af{i}((nx+np+1):end),nx+np,[]);
        end
    end
    for i = 1:nth
        if(second)
            ddrdu{i} = ddrdu{i}+fx_hessian_mex([xf;p0],nu,cases)'*bf{i}(1:(nx+np))*dtdu;
            bf{i}((nx+np+1):end) = reshape(ddrdu{i},[],1);
            if(length(tspan)>1&&all(-diff(tspan)>0))
                [~,boutk{i}] = ode45_rp(@ode1_a_hessian_mex,tspan,[bf{i};rtf{i}],1e-9,1e-12,[],tout,yout,deltat,nu,cases);
            else
                boutk{i} = [bf{i};rtf{i}];
            end
        else
            if(length(tspan)>1&&all(-diff(tspan)>0))
                [~,boutk{i}] = ode45_rp(@ode0_a_hessian_mex,tspan,[bf{i};rtf{i}],1e-9,1e-12,[],tout,yout,deltat,nu,cases);
            else
                boutk{i} = [bf{i};rtf{i}];
            end
        end
        bf{i} = boutk{i}(1:end-1,end);
        rtf{i} = boutk{i}(end,end);
        if(~fast)
            boutk{i} = boutk{i}(1:end-1,:);
            bout{i}(:,idx_tout) = fliplr(boutk{i});
        end
        if(second)
            ddrdu{i} = reshape(bf{i}((nx+np+1):end),nx+np,[]);
        end
    end
    tf = tspan(end);
    if(k>1)
        xf = xlim(:,k-1);
        if(iev(k-1))
            dhadxp = (iev(k-1)==(1:nh))*ghx_hessian_mex([xf;p0],nu,nh,civ(k-1));
            dhadt = dhadxp*f_hessian_mex([xf;p0],nu,civ(k-1));
        end
        if(second)
            if(iev(k-1))
                dtdu = -dhadxp*reshape(binterp_mex(tout,dxpduout,tf),nx+np,[])/dhadt;
            elseif(ismember(tf,ti))
                dtdu = zeros(1,nt+nx+np);
            else
                dtdu = [tv==tf,zeros(1,nx+np)];
            end
            for i = 1:nphi
                dcqdu{i} = dcqdu{i}-fx_hessian_mex([xf;p0],nu,cases)'*af{i}(1:(nx+np))*dtdu;
            end
            for i = 1:nth
                ddrdu{i} = ddrdu{i}-fx_hessian_mex([xf;p0],nu,cases)'*bf{i}(1:(nx+np))*dtdu;
            end
        end
        if(iev(k-1))
            fp = f_hessian_mex([xf;p0],nu,civ(k-1));
            fn = f_hessian_mex([xf;p0],nu,civ(k));
            if(second)
                dxpdu = [dxdulim(:,:,k-1);dpdu];
                d2hadtdu = dhadxp*fx_hessian_mex([xf;p0],nu,civ(k-1))*dxpdu;
                dfdxpp = fx_hessian_mex([xf;p0],nu,civ(k-1));
                dfdxpn = fx_hessian_mex([xf;p0],nu,civ(k));
            end
            for i = 1:nphi
                cq{i} = af{i}(1:(nx+np));
                dJdth{i} = (fp-fn)'*cq{i};
                af{i}(1:(nx+np)) = af{i}(1:(nx+np))-dhadxp'*dJdth{i}/dhadt;
                if(second)
                    d2Jdthdu{i} = (fp-fn)'*dcqdu{i}+cq{i}'*(dfdxpp-dfdxpn)*dxpdu;
                    dcqdu{i} = dcqdu{i}-dhadxp'*d2Jdthdu{i}/dhadt+(dhadxp'*dJdth{i}*d2hadtdu)/dhadt^2;
                end
            end
            for i = 1:nth
                dr{i} = bf{i}(1:(nx+np));
                dthdth{i} = (fp-fn)'*dr{i}+(th{i}==tf);
                bf{i}(1:(nx+np)) = bf{i}(1:(nx+np))-dhadxp'*dthdth{i}/dhadt;
                if(second)
                    d2thdthdu{i} = (fp-fn)'*ddrdu{i}+dr{i}'*(dfdxpp-dfdxpn)*dxpdu;
                    ddrdu{i} = ddrdu{i}-dhadxp'*d2thdthdu{i}/dhadt+(dhadxp'*dthdth{i}*d2hadtdu)/dhadt^2;
                end
            end
        end
        cases = civ(k-1);
    end
    for i = 1:nphi
        clim{i} = cat(2,af{i}(1:nx),clim{i});
        qlim{i} = cat(2,af{i}((nx+1):(nx+np)),qlim{i});
        if(second)
            af{i}((nx+np+1):end) = reshape(dcqdu{i},[],1);
            dcdulim{i} = cat(3,dcqdu{i}(1:nx,:),dcdulim{i});
            dqdulim{i} = cat(3,dcqdu{i}((nx+1):(nx+np),:),dqdulim{i});
        end
    end
    for i = 1:nth
        dlim{i} = cat(2,bf{i}(1:nx),dlim{i});
        rlim{i} = cat(2,bf{i}((nx+1):(nx+np)),rlim{i});
        if(second)
            bf{i}((nx+np+1):end) = reshape(ddrdu{i},[],1);
            dddulim{i} = cat(3,ddrdu{i}(1:nx,:),dddulim{i});
            drdulim{i} = cat(3,ddrdu{i}((nx+1):(nx+np),:),drdulim{i});
        end
    end
end
tlim = tlim(2:end);
nbt = length(tlim);
iev = or(iev,ismember(tlim,ti));
xplim = [xlim;repmat(p0,[1,nbt])];
xlim = xlim(:,~iev);
if(second)
    dxpdulim = [dxdulim;repmat(dpdu,[1,1,nbt])];
    dxdulim = dxdulim(:,:,~iev);
end
for i = 1:nphi
    cqout{i} = aout{i}(1:(nx+np),:);
    cqlim{i} = [clim{i};qlim{i}];
    clim{i} = clim{i}(:,[true,~iev(1:end)]);
    qlim{i} = qlim{i}(:,[true,~iev(1:end)]);
    if(second)
        dcqdulim{i} = [dcdulim{i};dqdulim{i}];
        dcdulim{i} = dcdulim{i}(:,:,[true,~iev(1:end)]);
        dqdulim{i} = dqdulim{i}(:,:,[true,~iev(1:end)]);
    end
end
for i = 1:nth
    drout{i} = bout{i}(1:(nx+np),:);
    drlim{i} = [dlim{i};rlim{i}];
    dlim{i} = dlim{i}(:,[true,~iev(1:end)]);
    rlim{i} = rlim{i}(:,[true,~iev(1:end)]);
    if(second)
        ddrdulim{i} = [dddulim{i};drdulim{i}];
        dddulim{i} = dddulim{i}(:,:,[true,~iev(1:end)]);
        drdulim{i} = drdulim{i}(:,:,[true,~iev(1:end)]);
    end
end
civp = civ(~iev);
civn = civ([false,~iev(1:end-1)]);
for i = 1:nth
    dthdts{i} = zeros(1,nt-1);
    for k = 1:(nt-1)
        xf = xlim(:,k);
        dthdts{i}(k) = [dlim{i}(:,k+1);rlim{i}(:,k+1)]'*(f_hessian_mex([xf;p0],nu,civp(k))-f_hessian_mex([xf;p0],nu,civn(k)));
    end
    dthdtf{i} = 0;
    dthdx0{i} = dlim{i}(:,1)';
    dthdp0{i} = rlim{i}(:,1)';
end
for i = 1:nphi
    dJdts{i} = zeros(1,nt-1);
    for k = 1:(nt-1)
        xf = xlim(:,k);
        dJdts{i}(k) = [clim{i}(:,k+1);qlim{i}(:,k+1)]'*(f_hessian_mex([xf;p0],nu,civp(k))-f_hessian_mex([xf;p0],nu,civn(k)))+lam(i,k+1);
    end
    xf = xlim(:,nt);
    dJdtf{i} = dJdxt{i}(1:(nx+1))*[1;eye(nx,nx+np)*f_hessian_mex([xf;p0],nu,civp(nt))]+lam(i,nt+1);
    dJdx0{i} = clim{i}(:,1)'+lam(i,nt+1+(1:nx));
    dJdp0{i} = qlim{i}(:,1)'+lam(i,nt+nx+1+(1:np));
end
if(second)
    for i = 1:nth
        d2thdtsdu{i} = zeros(nt-1,nt+nx+np);
        for k = 1:(nt-1)
            xf = xlim(:,k);
            d2thdtsdu{i}(k,:) = (f_hessian_mex([xf;p0],nu,civp(k))-f_hessian_mex([xf;p0],nu,civn(k)))'*[dddulim{i}(:,:,k+1);drdulim{i}(:,:,k+1)]+...
                [dlim{i}(:,k+1);rlim{i}(:,k+1)]'*(fx_hessian_mex([xf;p0],nu,civp(k))-fx_hessian_mex([xf;p0],nu,civn(k)))*[dxdulim(:,:,k);dpdu];
        end
        d2thdtfdu{i} = zeros(1,nt+nx+np);
        d2thdx0du{i} = dddulim{i}(:,:,1);
        d2thdp0du{i} = drdulim{i}(:,:,1);
    end
    for i = 1:nphi
        d2Jdtsdu{i} = zeros(nt-1,nt+nx+np);
        for k = 1:(nt-1)
            xf = xlim(:,k);
            d2Jdtsdu{i}(k,:) = (f_hessian_mex([xf;p0],nu,civp(k))-f_hessian_mex([xf;p0],nu,civn(k)))'*[dcdulim{i}(:,:,k+1);dqdulim{i}(:,:,k+1)]+...
                [clim{i}(:,k+1);qlim{i}(:,k+1)]'*(fx_hessian_mex([xf;p0],nu,civp(k))-fx_hessian_mex([xf;p0],nu,civn(k)))*[dxdulim(:,:,k);dpdu];
        end
        xf = xlim(:,nt);
        d2Jdtfdu{i} = f_hessian_mex([xf;p0],nu,civp(nt))'*[dcdulim{i}(:,:,nt+1);dqdulim{i}(:,:,nt+1)]+...
            dJdxt{i}(1:(nx+1))*[zeros(1,nt+nx+np);eye(nx,nx+np)*fx_hessian_mex([xf;p0],nu,civp(nt))*[dxdulim(:,:,nt);dpdu]]+...
            d2Jdxt2{i}(1,:)*[[zeros(1,nt-1),eye(1,1+nx+np)];dxdulim(:,:,nt);dxidu];
        d2Jdx0du{i} = dcdulim{i}(:,:,1);
        d2Jdp0du{i} = dqdulim{i}(:,:,1);
    end
end
tlim = {tlim};
xplim = {xplim};
if(second)
    dtdulim = {dtdulim};
    dxpdulim = {dxpdulim};
end
end