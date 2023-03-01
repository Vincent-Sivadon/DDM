function [Sd,Sp,bd,G,e] = GetDualStructures(nd,ne,L,E,F,S)

ni = nd-1;
Ld = 1/nd;
g = 1e6;

% Fields
f  = zeros(ne,1,nd);
f(end,1,end) = F;

% Stiffness Matrices
Ks = StiffnessMatricesDual(ne,nd,Ld,E,S);

% Dual terms
Sd = zeros(ni,ni);
Sp = zeros(ni,ni);
bd = zeros(ni,1);
e  = [];
G = [];
for i=1:nd
    % Isolates subpart of the stiffness matrix of subdomain i
    [Kii,Kib,Kbi,Kbb] = GetSubK(Ks,i,ne,nd);

    % Penalization
    if i==1
        Kii(1,1) = Kii(1,1) + g;
    end

    % Assembly Operator of subdomain i
    Abi = DualAssemblyOp(i,nd,ni);
    Ai  = PrimalAssemblyOp(i,nd,ni);

    % Primal Schur Complement (subdomain i)
    Spi = Kbb - Kbi* inv(Kii) * Kib;
    
    % Dual Schur Complement (subdomain i)
    Sdi = pinv(Spi);

    % Rigid Body Modes
    Rbi = null(Spi,'r');

    % Dual Operator : G
    tmp = Abi*Rbi;
    if(size(tmp))
        G = [G,tmp];
    end

    % Primal Right Hand Side (subdomain i)
    bids = Bids(i,ne,nd);
    iids = Iids(bids,ne);
    bpi = f(bids,:,i) - Kbi * inv(Kii) * f(iids,:,i);

    % Dual Right Hand Side (subdomain i)
    bdi = Sdi*bpi;

    % Dual Right Hand Side (2)
    e = [e;Rbi'*bpi];

    % Dual Operators
    Sd = Sd + Abi * Sdi * Abi';
    Sp = Sp + Ai  * Spi * Ai';
    bd = bd + Abi * bdi;
end

end % function
