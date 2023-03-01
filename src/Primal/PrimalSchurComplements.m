function Sps = PrimalSchurComplements(Ks,nd,ne)
    Sps = cell(nd,1);
    for i=1:nd
        [Kii,Kib,Kbi,Kbb] = GetSubK(Ks,i,ne,nd);
        Sps{i} = Kbb - Kbi * (Kii\Kib);
    end
end
