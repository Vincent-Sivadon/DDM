function A = PrimalAssemblyOp(i,nd)
    if i==1
        A = [1;0];
    elseif i==nd
        A = [0;1];
    else
        A = eye(2);
    end
end
