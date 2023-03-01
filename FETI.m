addpath(genpath("src"));

nd = 4;
ne = 3;
ni = nd-1;
L  = 1;
E  = 2e4;
F = 10;
S = 0.05;

[Sd,Sp,bd,G,e] = GetDualStructures(nd,ne,L,E,F,S);

% Preconditioner Optimal
wSd = (1.0/2.0) * eye(ni) * Sp;

tmp = G*inv(G'*G);
P = eye(ni) - tmp*G';
lbda = - tmp*e;
r   = P' * (-bd - Sd*lbda);     % Residual
z   = P*wSd*r;                  % Preconditioned Residual

m=100; count=0;
d = zeros(ni,1,m);
p = zeros(ni,1,m);
beta = zeros(m);
d(:,:,1) = z;

for i=1:m
    if (norm(r)<1e-3)
        break;
    end
    count = count + 1;

    p(:,:,i) = P' * Sd * d(:,:,i);

    alpha = dot(r,d(:,:,i)) / dot(d(:,:,i),p(:,:,i));

    lbda = lbda + alpha*d(:,:,i);

    r = r - alpha*p(:,:,i);

    z = P * wSd * r;

    for j=1:i
        beta(j) = - dot(z,p(:,:,j)) / dot(d(:,:,j),p(:,:,j));      
    end

    d(:,:,i+1) = z;
    for j=1:i
        d(:,:,i+1) = d(:,:,i+1) + beta(j)*d(:,:,i);
    end
end

alpha = inv(G'*G)*G'*(-bd - Sd*lbda);

% Isolates interior forces
Ab = [];
for i=1:nd
    Ab = [Ab,DualAssemblyOp(i,nd,ni)];
end
lbda_b = Ab'*lbda;

lbda_b
