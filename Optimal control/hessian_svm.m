function w = hessian_svm(Xs,Yt,p_mon)
m = size(Xs,1);
p_s = size(p_mon,1);
lambda = eps;
d = 1/m;
H12 = diag(Yt)*cell2mat(arrayfun(@(i)prod(Xs.^p_mon(i,:),2),1:p_s,'UniformOutput',false));
H = sparse([zeros(m,m),zeros(m,p_s);zeros(p_s,m),eye(p_s)]);
f = [-d*ones(m,1);zeros(p_s,1)];
A = [];
b = [];
Aeq = [H12',-eye(p_s)];
beq = zeros(p_s,1);
lb = [zeros(m,1);-Inf*ones(p_s,1)];
ub = [ones(m,1)/2/m/lambda;Inf*ones(p_s,1)];
x0 = [];
options = optimoptions('quadprog','Display','none');
x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options);
w = x((m+1):end);
idx_0 = ismember(p_mon,zeros(1,size(p_mon,2)),'rows');
w(idx_0) = w(idx_0)+d;
end
% max sum(z_i) / n + k w? w + sum(m_i (d - z_i - y_i w? f(x_i))) - sum(n_i z_i)
% 0 = 1 / n - m_i - n_i, i = 1,?,n
% 0 = 2 k w - sum(m_i y_i f(x_i))
% 0 <= m_i, i = 1,?,n
% 0 <= n_i, i = 1,?,n
% 
% max k w? w + sum(m_i d) - 2 k w? w
% 2 k w = sum(m_i y_i f(x_i))
% 0 <= m_i, i = 1,?,n
% m_i <= 1 / n, i = 1,?,n