addpath(genpath("src"))

nd = 10;      % Nb of subdomains
ne = 10;     % Nb of elements per subdomain
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
Abs = DualAssemblyOps(nd,ni);
Ks  = StiffnessMatrices(ne,nd,Ld,E,S);
Sps = PrimalSchurComplements(Ks,nd,ne);
bps = PrimalRHSs(Ks,f,nd,ne);

% Dual data structures
Sd = zeros(ni,ni);
bd = zeros(ni,1);
e = [];
G = [];

% Dual Operators
for i=1:nd
    [Kii,Kib,Kbi,Kbb] = GetSubK(Ks,i,ne,nd);

    % Rigib Body Modes
    Rbi = null(Sps{i},'r');

    % Dual Operator : G
    tmp = Abs{i}*Rbi;
    if (size(tmp))
        G = [G,tmp];
    end

    % Dual Schur Complement (subdomain i)
    Sdi = pinv(Sps{i});

    % Dual Right Hand Side (subdomain i)
    bdi = Sdi * bps{i};

    % Dual Right Hand Side (2)
    e = [e;Rbi'*bps{i}];

    % Dual Operators
    Sd = Sd + Abs{i} * Sdi * Abs{i}';
    bd = bd + Abs{i} * bdi;
end

% Assemble
Sgp = zeros(ni,ni);
bgp = zeros(ni,1);
for i=1:nd
    Sgp = Sgp + As{i} * Sps{i} * As{i}'; 
    bgp = bgp + As{i} * bps{i};
end

% Preconditionner
% wSp = (1/2)*eye(ni) * Sd * 2 * eye(ni);
wSp = Sd;

% Projector
tmp = G*inv(G'*Sgp*G);
P = eye(ni) - tmp*G'*Sgp;
ug = tmp*G'*bgp;
rg   = P' * bgp;               % Residual
z    = wSp * rg;               % Preconditioned Residual

m=10; count=0;
p = cell(m,1);
d = cell(m,1); d{1} = z;
beta = cell(m,1);
for i=1:m
    if (norm(rg)<1e-3)
        break;
    end
    count = count + 1;

    p{i} = Sp*P*d{i};

    % Optimal Step
    alpha = dot(rg,d{i}) / dot(d{i},p{i});

    % Iterate
    ug = ug + alpha*d{i};

    % Residual
    r = r - alpha*p{i};

    z = wSp * rg;

    % Reorthogonalization
    for j=1:i
        beta{j} = - dot(z,p{j}) / dot(d{j},p{j});
    end

    % Update search direction
    s = zeros(length(d{i}),1);
    for j=1:i
        s = s + beta{j}*d{i};
    end
    d{k+1} = z + s;
end

count
ug
