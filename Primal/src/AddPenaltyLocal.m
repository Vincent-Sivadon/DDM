function [A,b] = AddPenaltyLocal(A,b,D,id)
    % Let a linear system [A]{x} = {b} with CLs [D]{x} = {β}
    % The penalty method changes the system to :
    % [[A] + g[D]ᵗ[D]]{x} = {b} + g[D]ᵗ{β} with g big enough
    % Note : D is a dictionnary that associate interface indices with values
    % ex : D={1->0, 3->2} means x₁=0 and x₃=2
    % id is the subomain id

    g = 1e6;
    k = keys(D);
    for i=k'
        j = abs(i-id)+1;
        if j>1, break, end
        A(j,j) = A(j,j) + g;
        b(j)   = b(j) + g*D(j);
    end
end
