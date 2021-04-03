function ax = steps_1(sys_ol,sys_clr0,sys_cld0,sys_clw0,c_cl,n_cl,d_u,d_y,te)
n_y = length(c_cl)+length(n_cl);
n_u = length(d_u);
cols = 2;
rows = ceil((n_y+n_u)/cols);
ax = cell(n_y+n_u,3*n_u+n_y);
for j = 1:n_u
    figure;
    for i = 1:n_y
        k = find([c_cl,n_cl]==i,1);
        if(~isempty(k))
            ax{i,j} = subplot(rows,cols,i);stepplot(ax{i,j},sys_ol(k,j)*d_u(j),te);
        end
    end
    for i = 1:n_u
        ax{n_y+i,j} = subplot(rows,cols,n_y+i);stepplot(ax{n_y+i,j},sys_ol(n_y+i,j)*d_u(j),te);
    end
    set(gcf,'Units','normalized','OuterPosition',[0.2,0.2,0.6,0.6]);
end
for j = 1:n_u
    figure;
    for i = 1:n_y
        k = find([c_cl,n_cl]==i,1);
        if(~isempty(k))
            ax{i,n_u+j} = subplot(rows,cols,i);stepplot(ax{i,n_u+j},sys_clr0(k,j)*d_y(j),te);
        end
    end
    for i = 1:n_u
        ax{n_y+i,n_u+j} = subplot(rows,cols,n_y+i);stepplot(ax{n_y+i,n_u+j},sys_clr0(n_y+i,j)*d_y(j),te);
    end
    set(gcf,'Units','normalized','OuterPosition',[0.2,0.2,0.6,0.6]);
end
for j = 1:n_u
    figure;
    for i = 1:n_y
        k = find([c_cl,n_cl]==i,1);
        if(~isempty(k))
            ax{i,2*n_u+j} = subplot(rows,cols,i);stepplot(ax{i,2*n_u+j},sys_cld0(k,j)*d_u(j),te);
        end
    end
    for i = 1:n_u
        ax{n_y+i,2*n_u+j} = subplot(rows,cols,n_y+i);stepplot(ax{n_y+i,2*n_u+j},sys_cld0(n_y+i,j)*d_u(j),te);
    end
    set(gcf,'Units','normalized','OuterPosition',[0.2,0.2,0.6,0.6]);
end
for j = 1:n_y
    figure;
    for i = 1:n_y
        k = find([c_cl,n_cl]==i,1);
        if(~isempty(k))
            ax{i,3*n_u+j} = subplot(rows,cols,i);stepplot(ax{i,3*n_u+j},sys_clw0(k,j)*d_y(j)*tf(1,[6,1]),te);
        end
    end
    for i = 1:n_u
        ax{n_y+i,3*n_u+j} = subplot(rows,cols,n_y+i);stepplot(ax{n_y+i,3*n_u+j},sys_clw0(n_y+i,j)*d_y(j)*tf(1,[6,1]),te);
    end
    set(gcf,'Units','normalized','OuterPosition',[0.2,0.2,0.6,0.6]);
end
end
