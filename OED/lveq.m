function dXdt = lveq(t,X,Theta,dm,munf,r,cf) 

np = size(Theta,1);
dXdt = zeros(6,np);
X = reshape(X,1,6,np);

mu = -munf(permute(X,[1,3,2]),1)^2+r^2/(1+r^2)*munf(permute(X,[1,3,2]),2);

if(mu<=0)
    d = arc_eval(t,dm); %Given current time find the arc
else
    d = max(0,cf(permute(X,[1,3,2])));
end

for jj = 1:np

dxdt = zeros(6,1);
x = X(:,:,jj);
theta = Theta(jj,:);

dxdt(1) =  x(1) - (1+theta(1))*x(1)*x(2) - 0.4*x(1)*d;
dxdt(2) = -x(2) + (1+theta(2))*x(1)*x(2) - 0.2*x(2)*d;

r = [dxdt(1);dxdt(2)];

diff_states = length(r);

%Augment by sensitivity equations
rth = zeros(length(r),length(theta));
rth(1,1) = -x(1)*x(2); 
rth(1,2) = 0.0;
rth(2,1) = 0.0;
rth(2,2) =  x(1)*x(2);

rx = zeros(length(r),length(r));
rx(1,1) =  1 - (1+theta(1))*x(2) - 0.4*d;
rx(1,2) = -(1 + theta(1))*x(1);
rx(2,1) =  (1 + theta(2))*x(2);
rx(2,2) = -1 + (1+theta(2))*x(1) - 0.2*d;

seq = zeros(diff_states,length(theta));
seq(1,1) = x(3);
seq(1,2) = x(4);
seq(2,1) = x(5);
seq(2,2) = x(6);

counter = diff_states ; %initialize sensitivities counter by the number of 
%differential states 

for i = 1:diff_states
    for j = 1:length(theta)
    counter = counter+1;
    dxdt(counter) = rth(i,j)+(rx(i,1)*seq(1,j)+rx(i,2)*seq(2,j));
    end
end

dXdt(:,jj) = dxdt;

end

% if(mu>0)
%     keyboard
% end

dXdt = dXdt(:);

end

