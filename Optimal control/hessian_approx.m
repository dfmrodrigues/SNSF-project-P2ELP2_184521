function [N,p_mon,p_s,d,dus] = hessian_approx(n,nci,av,np,nt,idx_x0,m,ts0,tf0,x00,p00,deriv,du,sc,cheby,lb,ub)
nix = length(idx_x0);
np = nix*(deriv>0)+np;
N = nt+nix+np;
p_mon = calc_mon(n,N);
p_s = size(p_mon,1);
if(nt==length(ts0))
    u0 = [ts0,x00(idx_x0)'];
else
    u0 = [ts0,tf0,x00(idx_x0)'];
end
if(deriv)
    u0 = [u0,p00];
end
tb = [0,ts0];
te = [ts0,tf0];
e = @(i)[zeros(1,i-1),eye(1,N-i+1)];
tsi = @(i)e(i);
x0i = @(i)e(nt+i);
p0i = @(i)e(nt+nix+i);
d = zeros(p_s,length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*(length(lb)+length(ub)));
for k = 1:(length(ts0)+1)
    d(ismember(p_mon,zeros(1,N),'rows'),k) = te(k)-tb(k);
    if(deriv&&av(k)>0)
        d(ismember(p_mon,zeros(1,N),'rows'),k) = d(ismember(p_mon,zeros(1,N),'rows'),k)-(ub(1)-lb(1))/30/deriv;
    end
    if(k<=nt)
        d(ismember(p_mon,tsi(k),'rows'),k) = sc(k);
    end
    if(k>1)
        d(ismember(p_mon,tsi(k-1),'rows'),k) = -sc(k-1);
    end
end
for k = 1:nt
    d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*k-1) = 0.5;
    d(ismember(p_mon,tsi(k),'rows'),length(ts0)+1+2*k-1) = 1;
    d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*k) = 0.5;
    d(ismember(p_mon,tsi(k),'rows'),length(ts0)+1+2*k) = -1;
end
for j = 1:nix
    d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+2*j-1) = 0.5;
    d(ismember(p_mon,x0i(j),'rows'),length(ts0)+1+2*nt+2*j-1) = 1;
    d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+2*j) = 0.5;
    d(ismember(p_mon,x0i(j),'rows'),length(ts0)+1+2*nt+2*j) = -1;
end
if(deriv)
    for j = 1:nix
        d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+2*nix+2*j-1) = 0.5;
        d(ismember(p_mon,p0i(j),'rows'),length(ts0)+1+2*nt+2*nix+2*j-1) = 1;
        d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+2*nix+2*j) = 0.5;
        d(ismember(p_mon,p0i(j),'rows'),length(ts0)+1+2*nt+2*nix+2*j) = -1;
    end
end
if(~isempty(lb))
    for j = 1:nix
        d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j-1) = u0(nt+j)-lb(1);
        d(ismember(p_mon,x0i(j),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j-1) = sc(nt+j);
        d(:,length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j) = d(:,length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j-1);
        k = find(av==j);
        if(deriv)
            d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j) =...
                d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j)+u0(nt+nix+j)*(te(k)-tb(k));
            d(ismember(p_mon,p0i(j),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j) = sc(nt+nix+j)*(te(k)-tb(k));
            if(k<=nt)
                d(ismember(p_mon,tsi(k),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j) = sc(k)*u0(nt+nix+j);
                d(ismember(p_mon,p0i(j)+tsi(k),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j) = sc(nt+nix+j)*sc(k);
            end
            if(k>1)
                d(ismember(p_mon,tsi(k-1),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j) = -sc(k-1)*u0(nt+nix+j);
                d(ismember(p_mon,p0i(j)+tsi(k-1),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j) = -sc(nt+nix+j)*sc(k-1);
            end
        else
            d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j) =...
                d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j)+p00(j)*(te(k)-tb(k));
            if(k<=nt)
                d(ismember(p_mon,tsi(k),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j) = sc(k)*p00(j);
            end
            if(k>1)
                d(ismember(p_mon,tsi(k-1),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j) = -sc(k-1)*p00(j);
            end
        end
    end
end
if(~isempty(ub))
    for j = 1:nix
        d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j-1) = ub(1)-u0(nt+j);
        d(ismember(p_mon,x0i(j),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j-1) = -sc(nt+j);
        d(:,length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j) = d(:,length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j-1);
        k = find(av==j);
        if(deriv)
            d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j) =...
                d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j)-u0(nt+nix+j)*(te(k)-tb(k));
            d(ismember(p_mon,p0i(j),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j) = -sc(nt+nix+j)*(te(k)-tb(k));
            if(k<=nt)
                d(ismember(p_mon,tsi(k),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j) = -sc(k)*u0(nt+nix+j);
                d(ismember(p_mon,p0i(j)+tsi(k),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j) = -sc(nt+nix+j)*sc(k);
            end
            if(k>1)
                d(ismember(p_mon,tsi(k-1),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j) = sc(k-1)*u0(nt+nix+j);
                d(ismember(p_mon,p0i(j)+tsi(k-1),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j) = sc(nt+nix+j)*sc(k-1);
            end
        else
            d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j) =...
                d(ismember(p_mon,zeros(1,N),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j)-p00(j)*(te(k)-tb(k));
            if(k<=nt)
                d(ismember(p_mon,tsi(k),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j) = -sc(k)*p00(j);
            end
            if(k>1)
                d(ismember(p_mon,tsi(k-1),'rows'),length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j) = sc(k-1)*p00(j);
            end
        end
    end
end
rng(42,'twister');
validu = @(du)(prod(du.^p_mon,2)'*d);
if(cheby)
    nnn = ceil((ceil((m)^(1/N))+1)/2)*2+1;
    mn = nnn^(nt+nix);
    basenn = dec2base(0:(mn-1),nnn)-'0'+1;
    basenn(basenn>10) = basenn(basenn>10)-7;
    ptsn = sin(pi/2*(-1:(2/(nnn-1)):1))'*du(1,1:(nt+nix));
    dusn = cell2mat(cellfun(@(c)arrayfun(@(e)ptsn(c(e),e),1:(nt+nix)),num2cell(basenn,2),'UniformOutput',false));
    valid = zeros(mn,1);
    for l = 1:mn
        valid(l) = sum(validu([dusn(l,:),zeros(1,np)]./sc)<0)==0;
    end
    dusn = dusn(logical(valid),:);
    mn = size(dusn,1);
    if(np>0)
        nnp = floor((floor((m/mn)^(1/np))+1)/2)*2-1;
        mp = nnp^np;
        if(nnp>=2)
            basenp = dec2base(0:(mp-1),nnp)-'0'+1;
            basenp(basenp>10) = basenp(basenp>10)-7;
        else
            basenp = ones(1,np);
        end
        mt = mn*mp;
        dust = zeros(mt,N);
        for l = 1:mn
            dust(((l-1)*mp+1):(l*mp),1:(nt+nix)) = ones(mp,1)*dusn(l,:);
            lbp = zeros(1,np);
            ubp = zeros(1,np);
            for j = 1:np
                validulbj = @(du)(prod(du.^p_mon,2)'*d(:,length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+2*j));
                validuubj = @(du)(prod(du.^p_mon,2)'*d(:,length(ts0)+1+2*nt+(1+(deriv>0))*2*nix+(2*nix)*length(lb)+2*j));
                lbp(j) = max(-0.5*sc(nt+nix+j),min(0,fzero(@(p)validulbj([dusn(l,:),zeros(1,j-1),p,zeros(nix-j)]./sc),0)))+1e-9;
                ubp(j) = min(0.5*sc(nt+nix+j),max(0,fzero(@(p)validuubj([dusn(l,:),zeros(1,j-1),p,zeros(nix-j)]./sc),0)))-1e-9;
            end
            ptsp = (lbp+ubp)/2+sin(pi/2*(-1:(2/(nnp-1)):1))'*(ubp-lbp)/2;
            dust(((l-1)*mp+1):(l*mp),(nt+nix+1):N) =...
                cell2mat(cellfun(@(c)arrayfun(@(e)ptsp(c(e),e),1:np),num2cell(basenp,2),'UniformOutput',false));
        end
        valid = zeros(mt,1);
        for l = 1:mt
            valid(l) = sum(validu(dust(l,:)./sc)<-1e-9)==0;
        end
        dust = dust(logical(valid),:);
    else
        mt = mn;
        dust = dusn;
    end
else
    mt = 0;
    dust = [];
end
dus = [dust;zeros(m-mt,N)];
for l = (mt+1):m
    valid = false;
    while(~valid)
        dul = rand(1,N).*(du(2,:)+du(1,:))-du(1,:);
        valid = sum(validu(dul./sc)<0)==0;
    end
    dus(l,:) = dul;
end
for k = 1:(length(ts0)+1)
    if(deriv&&av(k)>0)
        d(ismember(p_mon,zeros(1,N),'rows'),k) = d(ismember(p_mon,zeros(1,N),'rows'),k)-(ub(1)-lb(1))/3/deriv;
    end
end
end