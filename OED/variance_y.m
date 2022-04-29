function m_p = variance_y(y,cases,theta,w,gv,Jv,sigma_z)
m_p0 = 0;
m_p1 = 0;
m_p2 = 0;
for jj = 1:size(theta,1)
    if(cases==1)
        m_p0 = m_p0+w(jj);
        m_p1 = m_p1+w(jj)*theta(jj,:)';
        m_p2 = m_p2+w(jj)*theta(jj,:)'*theta(jj,:);
    elseif(cases==2)
        m_p0 = m_p0+w(jj)*exp(-(y-gv(jj))^2/sigma_z^2/2)...
            /((2*pi)^(1/2)*(sigma_z^2)^(1/2));
        m_p1 = m_p1+w(jj)*exp(-(y-gv(jj))^2/sigma_z^2/2)...
            /((2*pi)^(1/2)*(sigma_z^2)^(1/2))*theta(jj,:)';
        m_p2 = m_p2+w(jj)*exp(-(y-gv(jj))^2/sigma_z^2/2)...
            /((2*pi)^(1/2)*(sigma_z^2)^(1/2))*theta(jj,:)'*theta(jj,:);
    elseif(cases==3)
        m_p0 = m_p0+w(jj)*exp(-((y-gv(jj))/Jv(jj))^2/sigma_z^2/2)...
            /((2*pi)^(1/2)*(sigma_z^2)^(1/2))/Jv(jj);
        m_p1 = m_p1+w(jj)*exp(-((y-gv(jj))/Jv(jj))^2/sigma_z^2/2)...
            /((2*pi)^(1/2)*(sigma_z^2)^(1/2))/Jv(jj)*theta(jj,:)';
        m_p2 = m_p2+w(jj)*exp(-((y-gv(jj))/Jv(jj))^2/sigma_z^2/2)...
            /((2*pi)^(1/2)*(sigma_z^2)^(1/2))/Jv(jj)*theta(jj,:)'*theta(jj,:);
    end
end
if(m_p0>1e-100)
    m_p = m_p2-(1/m_p0)*(m_p1*m_p1');
else
    m_p = 0*m_p2;
end
end