function [yJ,flag,tout,xpout,cqout,drout,uout] = hessian_eval(nci,phicase,av,nt,idx_x0,ti,ts0,tf0,x00,p00,deriv,lam,pi0,sc,du,hermite)
nx = length(x00);
nix = length(idx_x0);
np = length(p00)-(1-(deriv>0))*nix;
N = nt+nix+np;
t_mon = calc_mon(hermite,N);
t_s = size(t_mon,1);
m = size(du,1);
nphi = length(phicase);
nb = length(nci)-2*size(av,1);
J = cell(m,nphi);
th = cell(m,nb*nt);
dJdts = cell(m,nphi);
dJdtf = cell(m,nphi);
dJdx0 = cell(m,nphi);
dJdp0 = cell(m,nphi);
dthdts = cell(m,nb*nt);
dthdtf = cell(m,nb*nt);
dthdx0 = cell(m,nb*nt);
dthdp0 = cell(m,nb*nt);
if(hermite==2)
    d2Jdtsdu = cell(m,nphi);
    d2Jdtfdu = cell(m,nphi);
    d2Jdx0du = cell(m,nphi);
    d2Jdp0du = cell(m,nphi);
    d2thdtsdu = cell(m,nb*nt);
    d2thdtfdu = cell(m,nb*nt);
    d2thdx0du = cell(m,nb*nt);
    d2thdp0du = cell(m,nb*nt);
end
addpath('./hessian');
for l = 1:m
    fprintf('Evaluating point %d...\n',l);
    tsl = ts0+du(l,1:length(ts0));
    if(nt==length(ts0))
        tfl = tf0;
    else
        tfl = tf0+du(l,nt);
    end
    x0l = x00;
    x0l(idx_x0) = x0l(idx_x0)+du(l,nt+(1:nix))';
    if(deriv)
        p0l = p00'+du(l,nt+nix+(1:np))';
    else
        p0l = p00';
    end
    if(hermite==2)
        [J(l,:),thl,dJdts(l,:),dJdtf(l,:),dJdx0(l,:),dJdp0(l,:),dthdtsl,dthdtfl,dthdx0l,dthdp0l,...
            ~,~,~,~,tout,xpout,cqout,drout,uout,d2Jdtsdu(l,:),d2Jdtfdu(l,:),d2Jdx0du(l,:),d2Jdp0du(l,:),d2thdtsdul,d2thdtfdul,d2thdx0dul,d2thdp0dul] =...
            hessian_calc(nci,phicase,av,ti,0,tsl,tfl,x0l,p0l,lam,pi0);
    elseif(nargout>2)
        [J(l,:),thl,dJdts(l,:),dJdtf(l,:),dJdx0(l,:),dJdp0(l,:),dthdtsl,dthdtfl,dthdx0l,dthdp0l,...
            ~,~,~,~,tout,xpout,cqout,drout,uout] =...
            hessian_calc(nci,phicase,av,ti,0,tsl,tfl,x0l,p0l,lam,pi0);
    else
        [J(l,:),thl,dJdts(l,:),dJdtf(l,:),dJdx0(l,:),dJdp0(l,:),dthdtsl,dthdtfl,dthdx0l,dthdp0l] =...
            hessian_calc(nci,phicase,av,ti,0,tsl,tfl,x0l,p0l,lam,pi0);
    end
    if(nb>0)
        for k = 1:nt
            if(isempty(thl))
                th(l,k) = {Inf};
                dthdts(l,k) = {Inf*dJdts{l,1}};
                dthdtf(l,k) = {Inf*dJdtf{l,1}};
                dthdx0(l,k) = {Inf*dJdx0{l,1}};
                dthdp0(l,k) = {Inf*dJdp0{l,1}};
                if(hermite==2)
                    d2thdtsdu(l,k) = {Inf*d2Jdtsdu{l,1}};
                    d2thdtfdu(l,k) = {Inf*d2Jdtfdu{l,1}};
                    d2thdx0du(l,k) = {Inf*d2Jdx0du{l,1}};
                    d2thdp0du(l,k) = {Inf*d2Jdp0du{l,1}};
                end
            else
                if(k<nt)
                    th(l,k) = {(thl{1}-tsl(k))/sc(nt)};
                    dthdts(l,k) = {(dthdtsl{1}-((1:(nt-1))==k))/sc(nt)};
                    dthdtf(l,k) = {dthdtfl{1}/sc(nt)};
                else
                    th(l,k) = {(thl{1}-tfl)/sc(nt)};
                    dthdts(l,k) = {dthdtsl{1}/sc(nt)};
                    dthdtf(l,k) = {(dthdtfl{1}-1)/sc(nt)};
                end
                dthdx0(l,k) = {dthdx0l{1}/sc(nt)};
                dthdp0(l,k) = {dthdp0l{1}/sc(nt)};
                if(hermite==2)
                    d2thdtsdu(l,k) = {d2thdtsdul{1}/sc(nt)};
                    d2thdtfdu(l,k) = {d2thdtfdul{1}/sc(nt)};
                    d2thdx0du(l,k) = {d2thdx0dul{1}/sc(nt)};
                    d2thdp0du(l,k) = {d2thdp0dul{1}/sc(nt)};
                end
            end
        end
    end
end
rmpath('./hessian');
if(nb>0)
    valid = ~any(cell2mat(J(:,2:nphi))>0,2);
    if(~any(valid))
        signth = sortrows(sign(cell2mat(th)));
        signth = signth(ceil(m/10),:);
    else
        validth = cell2mat(th(valid,:));
        [~,iJ] = min(sum(validth,2));
        signth = sign(validth(iJ,:));
    end
    lastth = find(signth>=0,1,'last');
    firstth = find(signth<0,1,'first');
    flag = signth([lastth,firstth]).*cell2mat(th(:,[lastth,firstth]))>=0;
else
    flag = false(m,0);
end
yJ = cell(1,nphi);
for i = 1:nphi
    Ji = cell2mat(J(:,i));
    dJdtsi = cell2mat(dJdts(:,i));
    dJdtfi = cell2mat(dJdtf(:,i));
    dJdx0i = cell2mat(dJdx0(:,i));
    dJdp0i = cell2mat(dJdp0(:,i));
    yJ{i} = zeros(m,t_s);
    yJ{i}(:,1) = Ji(:);
    for k = 1:length(ts0)
        yJ{i}(:,k+1) = dJdtsi(:,k)*sc(k);
    end
    if(nt>length(ts0))
        yJ{i}(:,nt+1) = dJdtfi(:)*sc(nt);
    end
    for k = 1:nix
        yJ{i}(:,nt+k+1) = dJdx0i(:,idx_x0(k))*sc(nt+k);
    end
    if(deriv)
        for k = 1:np
            yJ{i}(:,nt+nix+k+1) = dJdp0i(:,k)*sc(nt+nix+k);
        end
    end
    if(hermite==2)
        d2Jdtsdui = cell2mat(cellfun(@(e)shiftdim(e,-1),d2Jdtsdu(:,i),'UniformOutput',false));
        d2Jdtfdui = cell2mat(cellfun(@(e)shiftdim(e,-1),d2Jdtfdu(:,i),'UniformOutput',false));
        d2Jdx0dui = cell2mat(cellfun(@(e)shiftdim(e,-1),d2Jdx0du(:,i),'UniformOutput',false));
        d2Jdp0dui = cell2mat(cellfun(@(e)shiftdim(e,-1),d2Jdp0du(:,i),'UniformOutput',false));
        for j = 1:length(ts0)
            for k = j:length(ts0)
                yJ{i}(:,j*(N+(1-j)/2)+k+1) = (d2Jdtsdui(:,j,k)+d2Jdtsdui(:,k,j))/2*sc(j)*sc(k);
            end
            if(nt<length(ts0))
                yJ{i}(:,j*(N+(1-j)/2)+nt+1) = (d2Jdtsdui(:,j,nt)+d2Jdtfdui(:,:,j))/2*sc(j)*sc(nt);
            end
            for k = 1:nix
                yJ{i}(:,j*(N+(1-j)/2)+nt+k+1) = (d2Jdtsdui(:,j,nt+idx_x0(k))+d2Jdx0dui(:,idx_x0(k),j))/2*sc(j)*sc(nt+k);
            end
            for k = 1:np
                yJ{i}(:,j*(N+(1-j)/2)+nt+nix+k+1) = (d2Jdtsdui(:,j,nt+nx+k)+d2Jdp0dui(:,k,j))/2*sc(j)*sc(nt+nix+k);
            end
        end
        if(nt<length(ts0))
            yJ{i}(:,nt*(N+(1-nt)/2)+nt+1) = d2Jdtfdui(:,:,nt)*sc(nt)^2;
            for k = 1:nix
                yJ{i}(:,nt*(N+(1-nt)/2)+nt+k+1) = (d2Jdtfdui(:,:,nt+idx_x0(k))+d2Jdx0dui(:,idx_x0(k),nt))/2*sc(nt)*sc(nt+k);
            end
            for k = 1:np
                yJ{i}(:,nt*(N+(1-nt)/2)+nt+nix+k+1) = (d2Jdtfdui(:,:,nt+nx+k)+d2Jdp0dui(:,k,nt))/2*sc(nt)*sc(nt+nix+k);
            end
        end
        for j = 1:nix
            for k = j:nix
                yJ{i}(:,(nt+j)*(N+(1-(nt+j))/2)+nt+k+1) = (d2Jdx0dui(:,idx_x0(j),nt+idx_x0(k))+d2Jdx0dui(:,idx_x0(k),nt+idx_x0(j)))/2*sc(nt+j)*sc(nt+k);
            end
            for k = 1:np
                yJ{i}(:,(nt+j)*(N+(1-(nt+j))/2)+nt+nix+k+1) = (d2Jdx0dui(:,idx_x0(j),nt+nx+k)+d2Jdp0dui(:,k,nt+idx_x0(j)))/2*sc(nt+j)*sc(nt+nix+k);
            end
        end
        for j = 1:np
            for k = j:np
                yJ{i}(:,(nt+nix+j)*(N+(1-(nt+nix+j))/2)+nt+nix+k+1) = (d2Jdp0dui(:,j,nt+nx+k)+d2Jdp0dui(:,k,nt+nx+j))/2*sc(nt+nix+j)*sc(nt+nix+k);
            end
        end
    end
end
end