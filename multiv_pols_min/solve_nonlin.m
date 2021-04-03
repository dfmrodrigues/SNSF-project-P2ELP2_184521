function sol = solve_nonlin(C,u_mon)
n = size(u_mon,2);
[G,d] = rref(C,1e-9);
r = length(d);
E = eye(n);
u_n = zeros(n,r);
for j = 1:r
    ref = find(ismemberrows_mex(u_mon,u_mon(d(j),:)+E));
    if(length(ref)~=n)
        sol = [];
        return;
    end
    u_n(:,j) = ref;
end
z = rand(n,1);
N = zeros(r,r,n);
for i = 1:n
    N(:,:,i) = G(1:r,u_n(i,:))'*z(i)/sum(z);
end
[Q,~] = schur(sum(N,3));
sol = zeros(r,n);
for i = 1:n
    sol(:,i) = sum(z)/z(i)*diag(Q'*N(:,:,i)*Q);
end
end