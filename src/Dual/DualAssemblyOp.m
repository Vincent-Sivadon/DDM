function Ab = DualAssemblyOp(i,nd,ni)
    if i==1
        Ab = zeros(ni,1);
        Ab(1) = 1;
    elseif i==nd
        Ab = zeros(ni,1);
        Ab(end) = -1;
    else
        Ab = zeros(ni,2);
        Ab(i-1:i,:) = [-1 0;0  1];
    end
end
