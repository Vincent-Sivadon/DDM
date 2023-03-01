addpath(genpath("src"))

nd = 3;      % Nb of subdomains
ne = 3;     % Nb of elements per subdomain
ni = nd-1;   % Nb of interfaces
L  = 1; Ld = L/nd; Le = Ld/(ne-1);  % Total, subdomain, and element length
E  = 2e4;    % Young's module
Fd = 10;     % Traction force
S  = 0.05;   % Surface area

% D is a dictionnary that associate indices with valyes
% It represents boundary conditions on the displacement field
D = dictionary();
D(1) = 0;  % associate first node to value 0

% u(:,:,i) is the field {u1,u2...,une}
u = zeros(ne,1,nd); ug = zeros(ni,1);
f = zeros(ne,1,nd); fg = zeros(ni,1);

% CG structures
r = zeros(ne,1,nd); rg = zeros(ni,1);
d = zeros(ne,1,nd); dg = zeros(ni,1);

% Data structure to store products Sp * db locally and globally
Spxdbs = cell(nd);
Spxbdg = zeros(ni,1);

% Boundary conditions on nodal forces
f(end)  = Fd;

% Get local operators
As  = PrimalAssemblyOps(nd,ni);
Ks  = StiffnessMatrices(ne,nd,Ld,E,S);
Sps = PrimalSchurComplements(Ks,nd,ne);
bps = PrimalRHSs(Ks,f,nd,ne);

% Compute local residual
for i=1:nd
    bids = Bids(i,ne,nd);
    r(bids,:,i) = bps{i} - Sps{i}*u(bids,:,i);
end

% Distribute local residual
for i=1:nd
    bids = Bids(i,ne,nd);
    rg = rg + As{i}*r(bids,:,i);
end
dg(:) = rg(:);

m = 100; count = 0;
for k=1:m
    % Distribute d
    for i=1:nd
        bids = Bids(i,ne,nd);
        d(bids,:,i) = As{i}' * dg;
    end

    % Compute local matrix-vector product
    for i=1:nd
        bids = Bids(i,ne,nd);
        Spxdbs{i} = Sps{i} * d(bids,:,i);
    end

    % Accumulate local Spxdb
    Spxdbg(:) = 0;
    for i=1:nd
        Spxdbg = Spxdbg + As{i} * Spxdbs{i};
    end

    % Compute optimal step
    alpha = dot(rg,dg) / dot(dg,Spxdbg);
    ug = ug + alpha*dg;
    u = u + alpha*d;
    rg = rg - alpha*Spxdbg;
    if (norm(rg)<1e-3)
        break;
    end
    beta = - dot(rg,Spxdbg) / dot(dg,Spxdbg);
    dg = rg + beta*dg;
end

ug
