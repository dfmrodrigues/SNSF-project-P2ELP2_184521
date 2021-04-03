function AJ = hessian_regress(sc,du,hermite,p_mon)
N = size(du,2);
p_s = size(p_mon,1);
t_mon = calc_mon(hermite,N);
t_s = size(t_mon,1);
m = size(du,1);
AJ = zeros(t_s*m,p_s);
for k = 1:t_s
    for j = 1:p_s
        AJ(((k-1)*m+1):(k*m),j) = prod(arrayfun(@(p,t)prod(p:-1:(p-t+1)),p_mon(j,:),t_mon(k,:)),2)*...
            prod((du./sc).^max(0,p_mon(j,:)-t_mon(k,:)),2);
    end
end
end