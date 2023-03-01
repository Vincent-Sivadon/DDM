function As = PrimalAssemblyOps(nd,ni)
    As = cell(nd,1);

    for i=1:nd
        As{i} = PrimalAssemblyOp(i,nd,ni);
    end
end
