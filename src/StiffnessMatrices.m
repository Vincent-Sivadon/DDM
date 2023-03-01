function Ks = StiffnessMatricesDual(ne,nd,Ld,E,S)
    % Stiffness Matrices
    Ks = zeros(ne,ne,nd);
    for i=1:nd
        % Reorder stiffness matrix
        K = StiffnessMatrix(ne,Ld,E,S);

        % Penalization
        if (i==1)
            K(1,1) = K(1,1) + 1e6;
        end

        bids = Bids(i,ne,nd);
        iids = Iids(bids,ne);
        rids = [iids,bids];    % Reordered indices
        K = K(:,rids);
        K = K(rids,:);

        Ks(:,:,i) = K;
    end
end
