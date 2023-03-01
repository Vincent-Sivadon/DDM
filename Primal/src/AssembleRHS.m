function bgp = AssembleRHS(bps,A,ni,nd)
    % Assemble primal right hand side from local ones
    % and local primal assembly operators

    bgp = zeros(ni,1);
    for i=1:nd
        bgp = bgp + A(:,:,i) * bps(:,:,i);
    end
end
