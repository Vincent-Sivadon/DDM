addpath(genpath("src"));

% Parameters
nd = 3;     % Nb of subdomains
ne = 3;     % Nb of elements
ni = nd-1;
L  = 1;     % Length
E  = 2e4;   % Young's Modulus
F = 10;     % Force
S = 0.05;   % Surface area

% Global structures
[Sd,Sp,bd,G,e] = GetDualStructures(nd,ne,L,E,F,S);

% Assemble operators for the direct solving method
LHS = [Sd, G;G',zeros(2,2)];
RHS = [-bd;-e];
v = LHS\RHS;
lambda_gb = v(1:end-2);
alpha_b   = v(end-1:end);

% Isolates interior forces
Ab = [];
for i=1:nd
    Ab = [Ab,DualAssemblyOp(i,nd)];
end
lambda_b = Ab'*lambda_gb;

