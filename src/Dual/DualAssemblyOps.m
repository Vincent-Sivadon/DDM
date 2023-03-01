function Abs = DualAssemblyOps(nd,ni)
    Abs = cell(nd,1);

    for i=1:nd
        Abs{i} = DualAssemblyOp(i,nd,ni);
    end
end
