function A = PrimalAssemblyOperators(ni,nd)
    % Primal Assembly Operators
    % Each of the A(:,:,i) belongs to a given process associated with a subdomain
    A = zeros(ni,2,nd);
    for i=1:nd
        A(:,:,i) = PrimalAssemblyOp(ni,nd,i);
    end
end
