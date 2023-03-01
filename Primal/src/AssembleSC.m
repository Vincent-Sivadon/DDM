function Sgb = AssembleSC(Sps,A,ni,nd)
    % Assemble primal schur complement from local ones
    % and from local primal assembly operators

    Sgb = zeros(ni,ni);
    for i=1:nd
        Sgb = Sgb + A(:,:,i)*Sps(:,:,i)*A(:,:,i)';
    end
end
