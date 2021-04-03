function [f,J,C,d,diff_pols_eval,diff_pols_deval,diff_basis_eval,diff_basis_deval] = quad_obj(x,mtrain,nparam,max_order,ords,pols,dpols,sc)
poly_num = size(ords,1); 
w = x(1:mtrain); 
points = x(mtrain+1:end);
theta = reshape(points,[mtrain,nparam]);
pols_eval = zeros(mtrain,nparam,max_order+1);
pols_deval = zeros(mtrain,nparam,max_order+1);
basis_eval = zeros(mtrain,poly_num);
basis_deval = zeros(mtrain,nparam,poly_num);

rhs = zeros(1,poly_num);
rhs(1,1) = 1.0;

for i = 1:mtrain
    [pols_eval(i,:,:),pols_deval(i,:,:)] = pols_mex(theta(i,:));
    [basis_eval(i,:),basis_deval(i,:,:)] = basis_mex(nparam,poly_num,ords',pols_eval(i,:,:),pols_deval(i,:,:));
end
if(nargout>4)
    diff_pols_eval = pols_eval;
    diff_pols_deval = pols_deval;
    diff_basis_eval = basis_eval;
    diff_basis_deval = basis_deval;
    for i = 1:mtrain
        for j = 1:poly_num
            diff_basis_eval(i,j) = diff_basis_eval(i,j)-prod(arrayfun(@(k)pols_eval(i,k,ords(j,k)+1),1:nparam));
            for l = 1:nparam
                diff_basis_deval(i,l,j) = diff_basis_deval(i,l,j)-pols_deval(i,l,ords(j,l)+1)*...
                    prod(arrayfun(@(k)pols_eval(i,k,ords(j,k)+1),[1:(l-1),(l+1):nparam]));
            end
        end
        for k = 1:nparam
            for n = 0:max_order
                diff_pols_eval(i,k,n+1) = diff_pols_eval(i,k,n+1)-pols{n+1,k}(theta(i,k));
                diff_pols_deval(i,k,n+1) = diff_pols_deval(i,k,n+1)-dpols{n+1,k}(theta(i,k));
            end
        end
    end
end

w_basis_deval = w.*basis_deval;
w_basis_deval = reshape(w_basis_deval,[],poly_num);

fv = w'*basis_eval-rhs;
Jv = [basis_eval;w_basis_deval]';
if(sc)
    f = fv*fv';
    J = 2*Jv'*fv';
else
    f = fv;
    J = Jv;
end
C = basis_eval';
d = rhs';
end