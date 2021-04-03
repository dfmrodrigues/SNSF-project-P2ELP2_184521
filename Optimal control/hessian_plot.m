function hessian_plot(nci,phicase,av,ti,t0,ts,tf,x0,p0,lam,pi0,nu_ineq)
[J,th,~,~,~,~,~,~,~,~,~,~,~,~,tout,xpout,cqout,drout,uout] =...
    hessian_calc(nci,phicase,av,ti,t0,ts,tf,x0,p0,lam,pi0);
nphi = length(J);
nth = length(th);
nu = size(uout,1);
nxp = size(xpout,1);
figure,plot(tout,uout,'LineWidth',1.5);
set(gca,'FontSize',20,'Position',[0.18,0.135,0.78,0.78]);
xlabel('$t$','Interpreter','LaTeX','FontSize',22);
if(nu==1)
    ylabel('$u(t)$','Interpreter','LaTeX','FontSize',22);
else
    ylabel('${\bf u}(t)$','Interpreter','LaTeX','FontSize',22);
    labelsu = cell(nu,1);
    for j = 1:nu
        labelsu{j} = ['$u_{',num2str(j),'}(t)$'];
    end
    legend(labelsu,'Location','Best','Interpreter','LaTeX','FontSize',22);
    legend('boxoff');
end
figure,plot(tout,xpout,'LineWidth',1.5);
set(gca,'FontSize',20,'Position',[0.18,0.135,0.78,0.78]);
xlabel('$t$','Interpreter','LaTeX','FontSize',22);
if(nxp==1)
    ylabel('$x(t)$','Interpreter','LaTeX','FontSize',22);
else
    ylabel('${\bf x}(t)$','Interpreter','LaTeX','FontSize',22);
    labelsx = cell(nxp,1);
    for j = 1:nxp
        labelsx{j} = ['$x_{',num2str(j),'}(t)$'];
    end
    legend(labelsx,'Location','Best','Interpreter','LaTeX','FontSize',22);
    legend('boxoff');
end
if(isempty(nu_ineq))
    for i = 1:nphi
        figure,plot(tout,cqout{i},'LineWidth',1.5);
        set(gca,'FontSize',20,'Position',[0.18,0.135,0.78,0.78]);
        xlabel('$t$','Interpreter','LaTeX','FontSize',22);
        ylabel(['{\boldmath$\lambda$}$^{\chi_{',num2str(i),'}}(t)$'],'Interpreter','LaTeX','FontSize',22);
        if(nxp==1)
            ylabel(['$\lambda^{\chi_{',num2str(i),'}}(t)$'],'Interpreter','LaTeX','FontSize',22);
        else
            ylabel(['{\boldmath$\lambda$}$^{\chi_{',num2str(i),'}}(t)$'],'Interpreter','LaTeX','FontSize',22);
            labelsx = cell(nxp,1);
            for j = 1:nxp
                labelsx{j} = ['$\lambda_{',num2str(j),'}^{\chi_{',num2str(i),'}}(t)$'];
            end
            legend(labelsx,'Location','Best','Interpreter','LaTeX','FontSize',22);
            legend('boxoff');
        end
    end
else
    for i = 2:nphi
        cqout{1} = cqout{1}+nu_ineq(i-1)*cqout{i};
    end
    figure,plot(tout,cqout{1},'LineWidth',1.5);
    set(gca,'FontSize',20,'Position',[0.18,0.135,0.78,0.78]);
    xlabel('$t$','Interpreter','LaTeX','FontSize',22);
    if(nxp==1)
        ylabel('$\lambda(t)$','Interpreter','LaTeX','FontSize',22);
    else
        ylabel('{\boldmath$\lambda$}$(t)$','Interpreter','LaTeX','FontSize',22);
        labelsx = cell(nxp,1);
        for j = 1:nxp
            labelsx{j} = ['$\lambda_{',num2str(j),'}(t)$'];
        end
        legend(labelsx,'Location','Best','Interpreter','LaTeX','FontSize',22);
        legend('boxoff');
    end
end
for i = 1:nth
    figure,plot(tout,drout{i},'LineWidth',1.5);
    set(gca,'FontSize',20,'Position',[0.18,0.135,0.78,0.78]);
    xlabel('$t$','Interpreter','LaTeX','FontSize',22);
    if(nxp==1)
        ylabel(['$\lambda^{\eta_{',num2str(i),'}}(t)$'],'Interpreter','LaTeX','FontSize',22);
    else
        ylabel(['{\boldmath$\lambda$}$^{\eta_{',num2str(i),'}}(t)$'],'Interpreter','LaTeX','FontSize',22);
        labelsx = cell(nxp,1);
        for j = 1:nxp
            labelsx{j} = ['$\lambda_{',num2str(j),'}^{\eta_{',num2str(i),'}}(t)$'];
        end
        legend(labelsx,'Location','Best','Interpreter','LaTeX','FontSize',22);
        legend('boxoff');
    end
end
end