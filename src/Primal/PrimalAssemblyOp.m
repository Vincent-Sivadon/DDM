function A = PrimalAssemblyOp(i,nd,ni)
    if i==1
        A = zeros(ni,1);
        A(1) = 1;
    elseif i==nd
        A = zeros(ni,1);
        A(end) = 1;
    else
        A = zeros(ni,2);
        A(i-1:i,:) = eye(2);
    end
end
