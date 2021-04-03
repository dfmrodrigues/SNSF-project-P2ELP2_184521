function [tau,p_sol,d_sol,deg,varargout] = multiv_pols_min(ab_mon,a,b,c,w,idx,rad,deg,varargin)
[~,n] = size(ab_mon);
n_p = size(b,2);
if(nargin>10)
    deg_max = varargin{3};
else
    deg_max = Inf;
end
if(nargin>8)
    Ik = varargin{1};
    Kj = varargin{2};
else
    Ik = 1:n;
    Kj = ones(1,size(c,2)/n_p);
end
if(nargout>4)
    solver = 1;
else
    solver = 0;
end
[p,n_k] = size(Ik);
n_a = size(idx,1);
Kj = [Kj,1:p];
n_c = size(c,2)/n_p+p;
c_deg = 1;
if(isempty(c))
    c = zeros(nchoosek(2*c_deg+n_k,n_k),0);
end
while(nchoosek(2*c_deg+n_k,n_k)<size(c,1))
    c_deg = c_deg+1;
end
ck_mon = calc_mon(2*c_deg,n_k);
ck_s = size(ck_mon,1);
c = mat2cell(c,nchoosek(2*c_deg+n_k,n_k),(n_c-p)*ones(1,n_p));
flagc = cell(size(c));
c_s = cell(size(c));
for m = 1:n_p
    c{m} = [c{m},zeros(ck_s,p)];
    c{m}(ismemberrows_mex(ck_mon,[zeros(1,n_k);2*c_deg*eye(n_k)]),end-p+(1:p)) = [rad^(2*c_deg);-ones(n_k,1)]*ones(1,p);
    c{m} = c{m}./max(abs(c{m}),[],1);
    c{m} = c{m}(:);
    flagc{m} = (c{m}~=0);
    c{m} = c{m}(flagc{m});
    c_s{m} = length(c{m});
end
while(deg<=deg_max)
    uniq = (deg==deg_max);
    u0k_mon = calc_mon(deg,n_k);
    u0k_s = size(u0k_mon,1);
    s0k_mon = calc_mon(2*deg,n_k);
    s0k_s = size(s0k_mon,1);
    s0_mon = sparse(s0k_s*p,n);
    for k = 1:p
        s0_mon((s0k_s*(k-1)+1):(s0k_s*k),Ik(k,:)) = s0k_mon;
    end
    s0_mon = full(s0_mon);
    [s0_mon,~,is0_mon] = uniquerows_mex(s0_mon,cumsum(s0_mon,2,'reverse'),Ik);
    s0_s = size(s0_mon,1);
    refs0 = cell(1,p);
    for k = 1:p
        refs0{k} = zeros(u0k_s,u0k_s);
    end
    for i = 1:u0k_s
        memb_i = ismemberrows_mex(s0k_mon,ones(u0k_s,1)*u0k_mon(i,:)+u0k_mon);
        for k = 1:p
            refs0{k}(i,:) = is0_mon(find(memb_i)+(k-1)*s0k_s);
        end
    end
    u1k_mon = calc_mon(deg-c_deg,n_k);
    u1k_s = size(u1k_mon,1);
    s1k_mon = calc_mon(2*(deg-c_deg),n_k);
    s1k_s = size(s1k_mon,1);
    s1_mon = sparse(s1k_s*p,n);
    for k = 1:p
        s1_mon((s1k_s*(k-1)+1):(s1k_s*k),Ik(k,:)) = s1k_mon;
    end
    s1_mon = full(s1_mon);
    [s1_mon,~,is1_mon] = uniquerows_mex(s1_mon,cumsum(s1_mon,2,'reverse'),Ik);
    s1_s = size(s1_mon,1);
    refs1 = cell(1,p);
    for k = 1:p
        refs1{k} = zeros(u1k_s,u1k_s);
    end
    for i = 1:u1k_s
        memb_i = ismemberrows_mex(s1k_mon,ones(u1k_s,1)*u1k_mon(i,:)+u1k_mon);
        for k = 1:p
            refs1{k}(i,:) = is1_mon(find(memb_i)+(k-1)*s1k_s);
        end
    end
    uck_mon = calc_mon(deg-c_deg,n_k);
    uck_s = size(uck_mon,1);
    sck_mon = calc_mon(2*(deg-c_deg),n_k);
    sck_s = size(sck_mon,1);
    refsc = cell(1,n_c);
    for j = 1:n_c
        refsc{j} = zeros(uck_s,uck_s);
    end
    for i = 1:uck_s
        memb_i = ismemberrows_mex(sck_mon,ones(uck_s,1)*uck_mon(i,:)+uck_mon);
        for j = 1:n_c
            refsc{j}(i,:) = find(memb_i)+(j-1)*sck_s;
        end
    end
    Q0_s = u0k_s*u0k_s*p;
    Qc_s = uck_s*uck_s*n_c;
    [~,pab_mon] = sortrows(cumsum(ab_mon,2,'reverse'));
    ab_mon = ab_mon(pab_mon,:);
    a = a(pab_mon,:);
    b = b(pab_mon,:);
    memb_k = ismemberrows_mex(s0_mon,ab_mon);
    s0_a = sparse(s0_s,n_a*n_p);
    s0_b = sparse(s0_s,n_p);
    s0_a(memb_k,:) = a;
    s0_b(memb_k,:) = b;
    s0_0 = speye(s0_s);
    R0 = sparse(s0_0(:,reshape(cell2mat(refs0),[],1)'));
    s0_ci = cell(size(c));
    s0_cj = cell(size(c));
    s0_cs = cell(size(c));
    s0_c = cell(size(c));
    Rc = cell(size(c));
    for m = 1:n_p
        s0_ci{m} = zeros(sck_s*c_s{m},1);
        s0_cj{m} = zeros(sck_s*c_s{m},1);
        s0_cs{m} = zeros(sck_s*c_s{m},1);
        for l = 1:sck_s
            memb_l = ismemberrows_mex(s0k_mon,ones(ck_s,1)*sck_mon(l,:)+ck_mon);
            ind_l = is0_mon(find(memb_l)+(Kj-1)*s0k_s)+(((1:n_c)-1)*sck_s+l-1)*s0_s;
            refs_l = ((l-1)*c_s{m}+1):(l*c_s{m});
            [s0_ci{m}(refs_l),s0_cj{m}(refs_l)] = ind2sub([s0_s,sck_s*n_c],ind_l(flagc{m}));
            s0_cs{m}(refs_l) = c{m};
        end
        s0_c{m} = sparse(s0_ci{m},s0_cj{m},s0_cs{m},s0_s,sck_s*n_c);
        Rc{m} = sparse(s0_c{m}(:,reshape(cell2mat(refsc),[],1)'));
    end
    q = length(w);
    fprintf('Minimizing with polynomial basis up to degree %d:\n',deg);
    s0_A = sparse(s0_s*n_p,q);
    s0_barA = sparse(s0_s*n_p,(Q0_s+Qc_s)*n_p);
    for m = 1:n_p
        refs_m = (s0_s*(m-1)+1):(s0_s*m);
        s0_A(refs_m,idx(:,m)) = s0_a(:,(n_a*(m-1)+1):(n_a*m));
        s0_barA(refs_m,(Q0_s*(m-1)+1):(Q0_s*m)) = R0;
        s0_barA(refs_m,(Q0_s*n_p+Qc_s*(m-1)+1):(Q0_s*n_p+Qc_s*m)) = Rc{m};
    end
    if(q==1)
        maxb = full(max(max(abs(s0_b(2:end,:)),[],1)));
        s0_b = s0_b/maxb;
    else
        maxb = 1;
    end
    if(q==1&&n_p==1&&s0_A(1)==1&&sum(s0_A(2:end)~=0)==0)
        q = 0;
    end
    if(solver>0)
        [Q0,z,tau,status,si] = sdp_sedumi(s0_A,s0_barA,s0_b,n_p,q,w,u0k_s,uck_s,Q0_s,Qc_s,p,n_c);
    else
        [Q0,z,tau,status,si] = sdp_mosek(s0_A,s0_barA,s0_b,n_p,q,w,u0k_s,uck_s,Q0_s,Qc_s,p,n_c);
    end
    tau = tau*maxb;
    y = reshape(z,s0_s,n_p);
    fprintf('Status returned by the solver: %s\n',status);
    fprintf('Optimal value of the cost function: %.4d\n',w*tau);
    d_sol = [];
    p_sol = [];
    for m = 1:n_p
        d_sol_m = [];
        Lr0 = u0k_s;
        Lr1 = u1k_s;
        for k = 1:p
            L0_mk = reshape(y(refs0{k},m),u0k_s,u0k_s);
            L1_mk = reshape(y(refs1{k},m),u1k_s,u1k_s);
            try
                [~,LS0,LV0] = svd(L0_mk);
                [~,LS1,~] = svd(L1_mk);
            catch
                tau = Inf;
                p_sol = zeros(1,n);
                d_sol = zeros(1,n);
                deg = 0;
                return;
            end
            Lr0 = min(Lr0,sum(sum(and(LS0>1e-6,LS0>LS0(1)/1000))));
            Lr1 = min(Lr1,sum(sum(and(LS1>1e-6,LS1>LS1(1)/1000))));
            if(uniq)
                Lr = 1;
            else
                Lr = Lr0;
            end
            if(Lr==0)
                break;
            end
            if(Lr>Lr1)
                d_sol_mk = [];
            else
                LC0 = LV0(:,1:Lr)';
                d_sol_mk = solve_nonlin(LC0,u0k_mon);
            end
            if(isempty(d_sol_mk))
                d_sol_m = [];
                break;
            end
            d_sol_m(1:Lr,Ik(k,:)) = d_sol_mk;
            d_sol_m = d_sol_m(1:Lr,:);
        end
        fprintf('Dual solution: moment matrices with rank %d and %d\n',Lr0,Lr1);
        if(Lr==0)
            continue;
        end
        if(isempty(d_sol_m))
            d_sol = [];
            break;
        end
        if(Lr0~=Lr)
            fprintf('Warning: dual solution may not be unique\n');
            deg = 0;
        end
        d_sol = [d_sol;d_sol_m];
        p_sol_m = [];
        Qr0 = 0;
        for k = 1:p
            Q0_mk = Q0(:,:,(m-1)*p+k);
            try
                [~,QS0,QV0] = svd(Q0_mk);
            catch
                tau = Inf;
                p_sol = zeros(1,n);
                d_sol = zeros(1,n);
                deg = 0;
                return;
            end
            Qr0 = max(Qr0,sum(sum(or(QS0>1e-8,QS0>1000*QS0(end)))));
            if(uniq)
                Qr = u0k_s-1;
            else
                Qr = Qr0;
            end
            if(Qr==u0k_s)
                break;
            end
            if(u0k_s-Qr>Lr1)
                p_sol_mk = [];
            else
                QC0 = QV0(:,(Qr+1):end)';
                p_sol_mk = solve_nonlin(QC0,u0k_mon);
            end
            if(isempty(p_sol_mk))
                p_sol_m = [];
                break;
            end
            p_sol_m(1:(u0k_s-Qr),Ik(k,:)) = p_sol_mk;
            p_sol_m = p_sol_m(1:(u0k_s-Qr),:);
        end
        fprintf('Primal solution: coefficient matrix with rank deficiency %d\n',u0k_s-Qr0);
        if(Qr==u0k_s)
            continue;
        end
        if(isempty(p_sol_m))
            p_sol = [];
            break;
        end
        if(Qr0~=Qr)
            fprintf('Warning: primal solution may not be unique\n');
            deg = 0;
        end
        p_sol = [p_sol;p_sol_m];
        if(Lr0>Lr1||u0k_s-Qr0>Lr1)
            fprintf('Warning: solution is not guaranteed to be global\n');
            deg = 0;
        end
    end
    if(~isempty(p_sol)&&~isempty(d_sol)&&(Lr==1||isinf(deg_max)))
        break;
    end
    deg = deg+1;
end
if(solver>0)
    varargout{1} = si;
end
end