function Sp = PrimalSchurComplement(i,K)
    % Return the local primal schur complement of subdomain i 
    
    % Isolates substructures of K
    [Kii,Kib,Kbi,Kbb] = GetSubK(K);

    Sp = Kbb - Kbi*(Kii\Kib);
end
