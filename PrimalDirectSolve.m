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
u = zeros(ne,1,nd);
f = zeros(ne,1,nd);

% CG structures
r = zeros(ne,1,nd); rg = zeros(ni,1);
d = zeros(ne,1,nd); dg = zeros(ni,1);

% Boundary conditions on nodal forces
f(end)  = Fd;

% Get local operators
As  = PrimalAssemblyOps(nd,ni);
Ks  = StiffnessMatrices(ne,nd,Ld,E,S);
Sps = PrimalSchurComplements(Ks,nd,ne);
bps = PrimalRHSs(Ks,f,nd,ne);

% Assemble
Sgp = zeros(ni,ni);
bgp = zeros(ni,1);
for i=1:nd
    Sgp = Sgp + As{i} * Sps{i} * As{i}'; 
    bgp = bgp + As{i} * bps{i};
end

% Solve Primal Interface Problem (direct solver)
% ----------------------------------------------
ug = Sgp \ bgp;
% ----------------------------------------------

% Distribute us
for i=1:nd
    bids = Bids(i,ne,nd);
    u(bids,:,i) = As{i}' * ug;
end

% Solve Dirichlet problem
for i=1:nd
    [Kii,Kib,Kbi,Kbb] = GetSubK(Ks,i,ne,nd);
    bids = Bids(i,ne,nd);
    iids = Iids(bids,ne);
    u(iids,:,i) = inv(Kii) * (f(iids,:,i) - Kib * u(bids,:,i));
end

ug
