function bps = PrimalRHSs(Ks,f,nd)
    % Primal Right Hand Side
    bps = zeros(2,1,nd);
    for i=1:nd
        fi = f([2:end-1],:,i);
        fb = f([1,end],:,i);
        [Kii,Kib,Kbi,Kbb] = GetSubK(Ks(:,:,i));
        bps(:,:,i) = fb - Kbi * (Kii\fi);
    end
end
