function s = newton_method(dg,g,s,p)
[n,v] = size(dg);
h = zeros(n,1);
j = zeros(n,v);
for i = 1:p
    for l = 1:n
        h(l) = g{l}(s);
        for m = 1:v
            j(l,m) = dg{l,m}(s);
        end
    end
    delta = -linsolve(j,h);
    s = s+delta*0.1;
    if(norm(delta)*1e9 < norm(s))
        break;
    end
end