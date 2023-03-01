function K = StiffnessMatrix(ne,Ld,E,S)
    % Return the local stiffness matrix reordered such that
    % K = {Kii Kib ; Kbi Kbb}
    % ne : nb of elements of the subdomain 
    % Ld : subdomain length 
    % E : Young's Module

    % Compute the operator B
    % For a subdomain i and an element j we have :
    % Bⱼⁱ = [dΦⱼⁱ(1)/dx dΦⱼⁱ(2)/dx]
    % with the shapes functions Φ (with linear assumption) :
    % Φⱼⁱ(1) = x2ⱼⁱ - x / x2ⱼⁱ - x1ⱼⁱ
    % Φⱼⁱ(2) = x - x1ⱼⁱ / x1ⱼⁱ - x2ⱼⁱ
    % such that
    % Bⱼⁱ = [-1 / x2ⱼⁱ - x1ⱼⁱ  ;  1 / x1ⱼⁱ - x2ⱼⁱ]
    % We then have Bⱼⁱ = [-1 1]/Le for all j and i with Le element length 
    B = zeros(ne,2);
    for iel=1:ne
        B(iel,:) = [-1 1]/(Ld/(ne-1));
    end

    % Compute the stiffness matrix for subdomain id
    % For a subdomain we have
    % K = ∫ Bⱼᵀ * E * Bⱼ dx (on the interval [x1ⱼ, x2ⱼ])
    % Since in 1D with linear assumption, Bⱼ is constant, we have the following
    K = zeros(ne,ne);
    for i=1:ne-1
        K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + S * B(i,:)' * E *  B(i,:) * (Ld/(ne-1));
    end

    % Reorder the matrix so that K = [Kii,Kib;Kbi,Kbb]
    rids = [2:ne-1,1,ne];
    K = K(:,rids);
    K = K(rids,:);
end
