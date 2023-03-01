function [Kii,Kib,Kbi,Kbb] = GetSubK(K)
    % Isolates subparts of the stiffness matrix K

    ne = length(K);

    % Range of interior and boundary indices
    irange = 1:ne-2;
    brange = ne-1:ne;

    Kii = K(irange,irange);
    Kib = K(irange,brange);
    Kbi = K(brange,irange);
    Kbb = K(brange,brange);
end
