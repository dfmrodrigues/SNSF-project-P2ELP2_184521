function res = calc_mon(deg,n)
s = nchoosek(deg+n,n);
root = [nchoosek(1:(deg+n),n),(deg+n+1)*ones(s,1)];
res = zeros(s,n);
for i = 1:s
    for j = 1:n
        res(end-i+1,j) = root(i,j+1)-root(i,j)-1;
    end
end