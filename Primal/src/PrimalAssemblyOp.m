function A = PrimalAssemblyOp(ni,nd,i)
    % Primal assembly opererator A for the subdomain i
    % Each line of A represents an interface
    % Columns represents indices of the boundary nodes of the subdomain i 
    % In 1D, there's 2 boundary nodes : 1b and 2b 
    % A(i,j) = 1 -> node jb is connected to the interface i 
    % A(i,j) = 0 -> node jb isn't
    A = zeros(ni,2);

    A(i:i+1,:) = [1 0;0 1];
end
