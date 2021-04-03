function out_seq = hessian_seq(in_seq,step,nt)
if(step==1)
    if(nt==2)
        out_seq = hessian_seq([],2,nt);
    else
        out_seq = [hessian_seq([],2,nt);hessian_seq(-1,2,nt);hessian_seq(-2,2,nt)];
    end
elseif(step==2)
    out_seq = [hessian_seq([in_seq,1,-1],3,nt);hessian_seq([in_seq,1,-2],3,nt)];
elseif(step==3)
    if(length(in_seq)==nt)
        out_seq = in_seq;
    elseif(length(in_seq)==nt-1)
        out_seq = hessian_seq(in_seq,4,nt);
    elseif(length(in_seq)==nt-2)
        out_seq = hessian_seq(in_seq,2,nt);
    else
        out_seq = [hessian_seq(in_seq,2,nt);hessian_seq(in_seq,4,nt)];
    end
elseif(step==4)
    if(in_seq(end)==-2)
        in_seq = [in_seq,-1];
    elseif(in_seq(end)==-1)
        in_seq = [in_seq,-2];
    end
    if(length(in_seq)==nt)
        out_seq = in_seq;
    else
        out_seq = hessian_seq(in_seq,2,nt);
    end
end
end