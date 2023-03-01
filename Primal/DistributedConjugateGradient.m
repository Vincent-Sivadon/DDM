addpath(genpath("src"))

nd = 100;      % Nb of subdomains
ne = 100;     % Nb of elements per subdomain
ni = nd+1;   % Nb of interfaces
L  = 1; Ld = L/nd; Le = Ld/(ne-1);  % Total, subdomain, and element length
E  = 2e4;    % Young's module
Fd = 10;     % Traction force
S  = 0.05;   % Surface area

% D is a dictionnary that associate indices with valyes
% It represents boundary conditions on the displacement field
D = dictionary();
D(1) = 0;  % associate first node to value 0

% u(:,:,i) is the field {u1,u2...,une}
[u,us] = NonOverlappingFields(ne,nd,ni);
[f,fs] = NonOverlappingFields(ne,nd,ni);

% CG structures
[r,rs] = NonOverlappingFields(ne,nd,ni);
[d,ds] = NonOverlappingFields(ne,nd,ni);

% Data structure to store products Sp db locally and globally
Spxdb = zeros(2,1,nd);
Spxdbs = zeros(ni,1);

% Boundary conditions on nodal forces
f(end)  = Fd;
fs(end) = Fd;

% Get local operators
A   = PrimalAssemblyOperators(ni,nd);
Ks  = StiffnessMatrices(ne,nd,Ld,E,S);
Sps = PrimalSchurComplements(Ks,nd);
bps = PrimalRHSs(Ks,f,nd);

% Add Penalty to substructures
for i=1:nd
    [Sps(:,:,i),bps(:,:,i)] = AddPenaltyLocal(Sps(:,:,i),bps(:,:,i),D,i);
end



% Distributed Conjugate Gradient
% ------------------------------

% Initialization
us(:) = 0;

% Compute local matrix-vector product
% Distribute us contribution to local u
for i=1:nd
    u([1,end],1,i) = A(:,:,i)' * us;
end

% Solve Dirichlet problem
for i=1:nd
    [Kii,Kib,Kbi,Kbb] = GetSubK(Ks(:,:,i));
    fi = f([2:end-1],:,i);
    ub = u([1,end],:,i);
    u([2:end-1],:,i) = inv(Kii) * (fi - Kib*ub);
end

% Compute local matrix-vector product
for i=1:nd
    r([1,end],:,i) = bps(:,:,i) - Sps(:,:,i)*u([1,end],:,i);
end

% Compute global residual
for i=1:nd
    rs(:) = rs(:) + A(:,:,i)*r([1,end],:,i);
end
ds(:) = rs(:);

m = 100; count = 0;
for k=1:m
    count = count + 1;

    % Compute local concatenated vector
    % Distribute ds
    for i=1:nd
        d([1,end],1,i) = A(:,:,i)' * ds;
    end
    
    % Solve local Dirichlet problem 
    for i=1:nd
        [Kii,Kib,Kbi,Kbb] = GetSubK(Ks(:,:,i));
        db = d([1,end],:,i);
        d([2:end-1],:,i) = -inv(Kii)*Kib*db;
    end

    % Compute local matrix-vector product 
    for i=1:nd
        db = d([1,end],:,i);
        Spxdb(:,:,i) = Sps(:,:,i) * db;
    end
    
    % Compute global matrix-vector product 
    Spxdbs(:) = 0;
    for i=1:nd
        Spxdbs(:) = Spxdbs(:) + A(:,:,i) * Spxdb(:,:,i);
    end

    % Compute optimal step 
    alpha = dot(rs,ds) / dot(ds,Spxdbs);
    us = us + alpha*ds;
    u = u + alpha * d;
    rs = rs - alpha*Spxdbs;
    if (norm(rs)<1e-3)
        break;
    end
    beta = - dot(rs,Spxdbs) / dot(ds,Spxdbs);
    ds = rs + beta*ds;
end

count
% us
