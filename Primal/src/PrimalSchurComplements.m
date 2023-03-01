function Sps = PrimalSchurComplements(Ks,nd)
    % Primal Schur Complements
    Sps = zeros(2,2,nd);
    for i=1:nd
        Sps(:,:,i) = PrimalSchurComplement(i,Ks(:,:,i));
    end
end
