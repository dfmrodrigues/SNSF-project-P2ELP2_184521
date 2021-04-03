function f = arc_eval(tcurr,dm)
%% Multistep heaviside function 
%Given switching times and corresponding values 
%find the value within a multistep f

tarr = dm(:,end); %last column has the switch times 

for i = 1:length(tarr)-1
    if tcurr >= tarr(i) && tcurr <= tarr(i+1)
        idx = i;
        f = dm(idx,1) + dm(idx,2)*(tcurr - tarr(i));
    elseif tcurr >= tarr(end)
        f = dm(end,1) + dm(end,2)*(tcurr - tarr(end));
    end
    
end 

end

