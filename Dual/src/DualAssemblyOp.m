function Ab = DualAssemblyOp(i,nd)
    if i==1
        Ab = [1;0];
    elseif i==nd
        Ab = [0;-1];
    else
        Ab = [-1 0;0 1];
    end
end
