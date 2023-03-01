function Ks = StiffnessMatrices(ne,nd,Ld,E,S)
    % Stiffness Matrices
    Ks = zeros(ne,ne,nd);
    for i=1:nd
        Ks(:,:,i) = StiffnessMatrix(ne,Ld,E,S);
    end
end
