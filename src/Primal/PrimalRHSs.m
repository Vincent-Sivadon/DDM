function bps = PrimalRHSs(Ks,f,nd,ne)
    bps = cell(nd,1);
    for i=1:nd
        [Kii,Kib,Kbi,Kbb] = GetSubK(Ks,i,ne,nd);
        bids = Bids(i,ne,nd);
        iids = Iids(bids,ne);
        bps{i} = f(bids,:,i) - Kbi * inv(Kii) * f(iids,:,i);
    end
end

