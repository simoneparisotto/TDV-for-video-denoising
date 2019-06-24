function K = Kmat(u)

[M,N,T] = size(u);

%% Derivatives operator (along i,j,k where k is the time)
D1 = spdiags([-ones(M,1) ones(M,1)],[0 1],M,M);
D2 = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N);
if T>1
    D3 = spdiags([-ones(T,1) ones(T,1)],[0 1],T,T);
else
    D3 = 1;
end
% Neumann boundary conditions
D1(M,:) = 0;
D2(N,:) = 0;
D3(T,:) = 0;
D1 = kron(speye(T),kron(speye(N),D1)); % i
D2 = kron(speye(T),kron(D2,speye(M))); % j
D3 = kron(D3,kron(speye(N),speye(M))); % k

K{1} = D1;
K{2} = D2;
K{3} = D3;

