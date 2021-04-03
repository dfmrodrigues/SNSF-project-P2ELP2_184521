function [tau,p_sol,d_sol,deg,tel,fvalv,uv,N,p_mon,p_s,d,dus,flag,cphi,As,ys,y,Av,yv,u0,du,sc] = hessian_solve(tf,x0,nci,nphi,av,np,ni,nt,ti,scv,deriv,lb,ub,varargin)
if(length(av)==1)
    av = hessian_seq([],1,av);
end
nseq = size(av,1);
hermite = 1;
if(isempty(varargin))
    n = 6;
    m = 1000;
    cheby = false;
    suppl = 0;
else
    n = varargin{1};
    m = varargin{2};
    cheby = varargin{3};
    suppl = varargin{4};
    if(length(varargin)>4)
        lam = varargin{5};
        pi0 = varargin{6};
    else
        lam = zeros(nphi,size(av,2)+length(x0)+2*ni+1);
        pi0 = zeros(size(av,2)+length(x0)+2*ni,1);
    end
end
du = zeros(2,nt+ni);
sc = zeros(1,nt+ni);
ts0 = zeros(1,size(av,2)-1);
for k = 1:length(ts0)
    ts0(:,k) = tf/2+(k-(length(ts0)+1)/2)*tf/200;
    du(:,k) = [tf/2-((length(ts0)+1)/2-1)*tf/200;tf/2-((length(ts0)+1)/2-1)*tf/200];
    sc(k) = tf-((length(ts0)+1)/2-1)*tf/100;
end
if(nt==length(ts0))
    tf0 = tf;
else
    tf0 = 5*tf/6;
    du(:,nt) = [tf/6;tf/6];
    sc(nt) = tf/3;
end
x00 = [x0;(lb+ub)/2*ones(ni,1)];
du(:,nt+(1:ni)) = [(ub-lb)/2*ones(1,ni);(ub-lb)/2*ones(1,ni)];
sc(nt+(1:ni)) = (ub-lb)*ones(1,ni);
idx_x0 = length(x00)+((1-ni):0);
if(nt==length(ts0))
    u0 = [ts0,x00(idx_x0)'];
else
    u0 = [ts0,tf0,x00(idx_x0)'];
end
if(deriv)
    u0 = [u0,zeros(1,ni)];
    du = [du,[deriv*ones(1,ni);deriv*ones(1,ni)]];
    sc = [sc,2*deriv*ones(1,ni)];
end
tau = cell(nseq,1);
p_sol = cell(nseq,1);
d_sol = cell(nseq,1);
deg = cell(nseq,1);
tel = cell(nseq,1);
fvalv = cell(nseq,1);
uv = cell(nseq,1);
N = cell(nseq,1);
p_mon = cell(nseq,1);
p_s = cell(nseq,1);
d = cell(nseq,1);
dus = cell(nseq,1);
flag = cell(nseq,1);
cphi = cell(nseq,1);
As = cell(nseq,1);
ys = cell(nseq,1);
y = cell(nseq,1);
Av = cell(nseq,1);
yv = cell(nseq,1);
tout = cell(nseq,1);
xpout = cell(nseq,1);
cqout = cell(nseq,1);
drout = cell(nseq,1);
uout = cell(nseq,1);
delete(gcp('nocreate'));
if(nseq==1)
    workers = 0;
else
    workers = Inf;
end
parfor(l = 1:nseq,workers)
% for l = 1:nseq
maxNumCompThreads(2);
av_l = av(l,:);
p00 = zeros(1,ni);
if(~deriv)
    [u_l,fval_l] = hessian_fmincon(nci,nphi,av_l,np,idx_x0,ti,ts0,tf0,x00,p00,2*scv,lam,pi0,lb,ub,0.1);
    p00 = u_l((length(ts0)+1+ni+1):end);
    fprintf('Local optimization:\n');
    disp(fval_l);
    disp(u_l);
%     p00 = 0;
end
tstart_l = tic;
ts0_l = ts0;
tf0_l = tf0;
x00_l = x00;
p00_l = p00;
u0_l = u0;
du_l = du;
sc_l = sc;
[N{l},p_mon{l},p_s{l},d{l},dus{l}] = hessian_approx(n,nci,av_l,np,nt,idx_x0,m,ts0_l,tf0_l,x00_l,p00_l,deriv,du_l,sc_l,cheby,lb,ub);
fprintf('Function 1 to %d, sample points:\n',nphi);
[ys{l},flag{l}] = hessian_eval(nci,1:nphi,av_l,nt,idx_x0,ti,ts0_l,tf0_l,x00_l,p00_l,deriv,lam,pi0,sc_l,dus{l},hermite);
if(suppl)
    sys_l = sort(ys{l}{1}(:,1));
    sel_l = ys{l}{1}(:,1)<=sys_l(suppl);
    min_l = max(u0_l+min(dus{l}(sel_l,:)),u0_l-du_l(1,:));
    max_l = min(u0_l+max(dus{l}(sel_l,:)),u0_l+du_l(2,:));
    u0_l = (min_l+max_l)/2;
    du_l = [1;1]*(max_l-min_l)/2;
    sc_l = max(1e-3,max_l-min_l);
    ts0_l = u0_l(1:length(ts0_l));
    if(nt~=length(ts0_l))
        tf0_l = u0_l(nt);
    end
    x00_l(idx_x0) = u0_l(nt+(1:ni));
    if(deriv)
        p00_l = u0_l((nt+ni+1):end);
    end
    [N{l},p_mon{l},p_s{l},d{l},dus{l}] = hessian_approx(n,nci,av_l,np,nt,idx_x0,m,ts0_l,tf0_l,x00_l,p00_l,deriv,du_l,sc_l,cheby,lb,ub);
    fprintf('Function 1 to %d, supplementary sample points:\n',nphi);
    [ys{l},flag{l}] = hessian_eval(nci,1:nphi,av_l,nt,idx_x0,ti,ts0_l,tf0_l,x00_l,p00_l,deriv,lam,pi0,sc_l,dus{l},hermite);
end
for i = 1:size(flag{l},2)
    fprintf('Constructing support vector machine %d...\n',i);
    w_i = hessian_svm(dus{l}./sc_l,2*flag{l}(:,i)-1,p_mon{l});
    d{l} = [d{l},w_i];
end
flag{l} = all(flag{l},2);
% dus{l} = dus{l}(flag{l},:);
% for i = 1:nphi
%     ys{l}{i} = reshape(ys{l}{i}(flag{l},:),[],1);
% end
for i = 1:nphi
    ys{l}{i} = reshape(ys{l}{i},[],1);
end
data_s = matfile(['data_',num2str(av_l),'.mat'],'Writable',true);
data_s.dus = dus{l};
data_s.ys = ys{l};
data_s.flag = flag{l};
% data_s = matfile(['data_',num2str(av_l),'.mat'],'Writable',false);
% dus{l} = data_s.dus;
% ys{l} = data_s.ys;
% flag{l} = data_s.flag;
for i = 1:nphi
    ys{l}{i} = reshape(ys{l}{i},size(dus{l},1),[]);
end
As{l} = hessian_regress(sc_l,dus{l}(flag{l},:),hermite,p_mon{l});
% As{l} = As{l}(1:size(dus{l}(flag{l},:),1),:);
% M_l = eye(size(As{l},1))-As{l}*pinv(As{l});
cphi{l} = cell(1,nphi);
for i = 1:nphi
%     ys{l}{i} = ys{l}{i}(:,1);
    cphi{l}{i} = linsolve(As{l},reshape(ys{l}{i}(flag{l},:),[],1));
%     var_y = sum((reshape(ys{l}{i},m,[])-mean(reshape(ys{l}{i},m,[]),1)).^2,1);
%     w_y = reshape(repmat(1./sqrt(var_y),m,1),[],1);
%     cphi{l}{i} = linsolve(w_y.*As{l},w_y.*ys{l}{i});
%     e_loo_l_i = sum(reshape(((1./diag(M_l)).*M_l*ys{l}{i}).^2,m,[]),1)./var_y;
%     disp(e_loo_l_i)
end
cphi_s = matfile(['cphi_',num2str(av_l),'.mat'],'Writable',true);
cphi_s.cphi = cphi{l};
% cphi_s = matfile(['cphi_',num2str(av_l),'.mat'],'Writable',false);
% cphi{l} = cphi_s.cphi;
[tau{l},p_sol{l},d_sol{l},deg{l}] = hessian_optim(n,N{l},p_mon{l},p_s{l},d{l},cphi{l},u0_l,sc_l);
tsv_l = d_sol{l}(1:length(ts0_l));
if(nt==length(ts0_l))
    tfv_l = tf0_l;
else
    tfv_l = d_sol{l}(nt);
end
x0v_l = x00_l;
x0v_l(idx_x0) = d_sol{l}(nt+(1:ni));
if(~deriv)
    p_sol{l} = [p_sol{l},p00_l];
    d_sol{l} = [d_sol{l},p00_l];
end
p0v_l = d_sol{l}((nt+ni+1):end);
% uv{l} = d_sol{l};
[uv{l},fvalv{l}] = hessian_fmincon(nci,nphi,av_l,np,idx_x0,ti,tsv_l+1e-6,tfv_l+2e-6,x0v_l,p0v_l,2*scv,lam,pi0,lb,ub,0.001);
tel{l} = toc(tstart_l);
if(~deriv)
    duv_l = uv{l}([1:nt,(length(ts0_l)+2):(length(ts0_l)+1+ni)])-u0_l;
else
    duv_l = uv{l}([1:nt,(length(ts0_l)+2):end])-u0_l;
end
fprintf('Function 1 to %d, validation points:\n',nphi);
[yv{l},~,tout{l},xpout{l},cqout{l},drout{l},uout{l}] = hessian_eval(nci,1:nphi,av_l,nt,idx_x0,ti,ts0_l,tf0_l,x00_l,p00_l,deriv,lam,pi0,sc_l,duv_l,hermite);
Av{l} = hessian_regress(sc_l,duv_l,hermite,p_mon{l});
% fvalv{l} = yv{l}{1}(1);
y{l} = cell(1,length(cphi{l}));
for i = 1:length(cphi{l})
    y{l}{i} = Av{l}*cphi{l}{i};
    disp([y{l}{i},yv{l}{i}(:)]);
end
% teb = [0,uv{l}(1:(length(ts0_l)+1))];
% prev = [];
% figure;
% for k = 1:(length(ts0_l)+1)
%     if(av_l(k)>0)
%         j = av_l(k);
%         plot([teb(k),teb(k)],[prev,uv{l}((length(ts0_l)+1)+j)],'LineWidth',1.5);
%         hold on;ax = gca;ax.ColorOrderIndex = 1;
%         prev = uv{l}((length(ts0_l)+1)+j)+uv{l}((length(ts0_l)+1)+ni+j)*(teb(k+1)-teb(k));
%         plot([teb(k),teb(k+1)],[uv{l}((length(ts0_l)+1)+j),prev],'LineWidth',1.5);
%         hold on;ax = gca;ax.ColorOrderIndex = 1;
%     elseif(av_l(k)==-1)
%         plot([teb(k),teb(k)],[prev,lb(1)],'LineWidth',1.5);
%         hold on;ax = gca;ax.ColorOrderIndex = 1;
%         prev = lb(1);
%         plot([teb(k),teb(k+1)],[lb(1),prev],'LineWidth',1.5);
%         hold on;ax = gca;ax.ColorOrderIndex = 1;
%     elseif(av_l(k)==-2)
%         plot([teb(k),teb(k)],[prev,ub(1)],'LineWidth',1.5);
%         hold on;ax = gca;ax.ColorOrderIndex = 1;
%         prev = ub(1);
%         plot([teb(k),teb(k+1)],[ub(1),prev],'LineWidth',1.5);
%         hold on;ax = gca;ax.ColorOrderIndex = 1;
%     end
% end
% set(gca,'FontSize',20,'Position',[0.18,0.135,0.78,0.78]);
% xlabel('$t$ [min]','Interpreter','LaTeX','FontSize',22);
% ylabel('$u^{*}(t)$ [mL min$^{-1}$]','Interpreter','LaTeX','FontSize',22);
fprintf('Global optimization:\n');
disp(tau{l});
disp(p_sol{l});
disp(d_sol{l});
disp(deg{l});
fprintf('Global+local optimization:\n');
disp(tel{l});
disp(fvalv{l});
disp(uv{l});
end
delete(gcp('nocreate'));
lmin = min(cell2mat(fvalv))==cell2mat(fvalv);
figure,plot(tout{lmin},uout{lmin},'LineWidth',1.5);hold on;
set(gca,'FontSize',20,'Position',[0.18,0.135,0.78,0.78]);
xlabel('$t$ [min]','Interpreter','LaTeX','FontSize',22);
ylabel('$u^{*}(t)$ [mL min$^{-1}$]','Interpreter','LaTeX','FontSize',22);
end