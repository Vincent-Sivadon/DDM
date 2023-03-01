addpath(genpath("src"))

nd = 80;    % Nb of subdomains
ne = 800;    % Nb of elements per subdomain
ni = nd+1;  % Nb of interfaces
L  = 1; Ld = L/nd; Le = Ld/(ne-1);  % Total, subdomain, and element lengths
E  = 2e4;   % Young's module
Fd = 10;    % Traction force
S  = 0.05;  % Surface area

% D is a dictionnary that associate indice with values
% It represents boundary conditions on the displacement field
D = dictionary();
D(1) = 0;  % associate first node to value 0

% u(:,:,i) is the field {u1,u2...,une}
[u,us] = NonOverlappingFields(ne,nd,ni);
[f,fs] = NonOverlappingFields(ne,nd,ni);

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

% Get global interface operators
Sgb = AssembleSC(Sps,A,ni,nd);
bgp = AssembleRHS(bps,A,ni,nd);

% Contionning informations
c = cond(Sgb)

% Direct Solving
% -----------------
tic
    us = Sgb \ bgp;
toc
% -----------------

% Distribute us
for i=1:nd
    u([1,end],:,i) = A(:,:,i)' * us;
end

% Solve Dirichlet problem
for i=1:nd
    [Kii,Kib,Kbi,Kbb] = GetSubK(Ks(:,:,i));
    u([2:end-1],:,i) = inv(Kii) * (f([2:end-1],:,i) - Kib * u([1,end],:,i));
end
 
% us
