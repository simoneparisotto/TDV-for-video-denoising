function AA = compute_average_order(M,N,T)

%% B is the matrix for averaging to bring the derivative in the cell centres

%% FOR ORDER 1
% (u_{i,j,k} + u_{i+1,j,k})/2
B1plus                  = 0.5*spdiags([ones(M,1) ones(M,1)],[0 1],M,M);
B1plus(end,[end-1 end]) = [0.5 0.5];
B{1}                    = kron(speye(T),kron(speye(N),B1plus));

% (u_{i,j,k} + u_{i,j+1,k})/2
B2plus                  = 0.5*spdiags([ones(N,1) ones(N,1)],[0 1],N,N);
B2plus(end,[end-1 end]) = [0.5 0.5];
B{2}                    = kron(speye(T),kron(B2plus,speye(M)));

% (u_{i,j,k} + u_{i,j,k+1})/2
B3plus                  = 0.5*spdiags([ones(T,1) ones(T,1)],[0 1],T,T);
B3plus(end,[end-1 end]) = [0.5 0.5];
B{3}                    = kron(B3plus,kron(speye(N),speye(M)));

%% FOR ORDER 2
% (u_{i,j,k} + u_{i-1,j+1,k})/2
B4plus                  = 0.5*spdiags([ones(M*N,1) ones(M*N,1)],[0 M-1],M*N,M*N);
B{4}                    = kron(speye(T),B4plus);
% (u_{i,j,k} + u_{i+1,j-1,k})/2
B5plus                  = 0.5*spdiags([ones(M*N,1) ones(M*N,1)],[0 -M+1],M*N,M*N);
B{5}                    = kron(speye(T),B5plus);
%(u_{i,j,k} + u_{i-1,j,k+1})/2
B6plus                  = 0.5*spdiags([ones(M*N*T,1) ones(M*N*T,1)],[0 M*N-1],M*N*T,M*N*T);
B{6}                    = B6plus;
% (u_{i,j,k} + u_{i+1,j,k-1})/2
B7plus                  = 0.5*spdiags([ones(M*N*T,1) ones(M*N*T,1)],[0 -M*N+1],M*N*T,M*N*T);
B{7}                    = B7plus;
% (u_{i,j,k} + u_{i,j-1,k+1})/2
B8plus                  = 0.5*spdiags([ones(M*N*T,1) ones(M*N*T,1)],[0 (N-1)*T],M*N*T,M*N*T);
B{8}                    = B8plus;
% (u_{i,j,k} + u_{i,j+1,k-1})/2
B9plus                  = 0.5*spdiags([ones(M*N*T,1) ones(M*N*T,1)],[0 -(N-1)*T],M*N*T,M*N*T);
B{9}                    = B9plus;

%% FOR ORDER 3


%% R,L we need to delete rows and columns; R for rows L for columns

% I-J (2:end-2,2:end-2,2:end-2)
r5                    = ones(M,N,T); 
r5([1,end-1,end],:,:) = 0; 
r5                    = find(r5(:));
l5                    = ones(M-3,N,T); 
l5(:,[1,end-1,end],:) = 0; 
l5                    = find(l5(:));
k5                    = ones(M-3,N-3,T); 
k5(:,:,[1,end-1,end]) = 0;
k5                    = find(k5(:));

% % I-K (2:end-2,2:end-2,2:end-2)
% r6                    = ones(M,N,T); 
% r6(:,[1,end-1,end],:) = 0; 
% r6                    = find(r6(:));
% l6                    = ones(M,N-3,T); 
% l6(:,:,[1,end-1,end]) = 0; 
% l6                    = find(l6(:));
% k6                    = ones(M,N-3,T-3); 
% k6([1,end-1,end],:,:) = 0;
% k6                    = find(k6(:));
% 
% % J-K (2:end-2,2:end-2)
% r7                    = ones(M,N,T); 
% r7(:,:,[1,end-1,end]) = 0; 
% r7                    = find(r7(:));
% l7                    = ones(M,N,T-3); 
% l7([1,end-1,end],:,:) = 0; 
% l7                    = find(l7(:));
% k7                    = ones(M-3,N,T-3); 
% k7(:,[1,end-1,end],:) = 0;
% k7                    = find(k7(:));
% 

% I-J
R{5} = spdiags(ones(M*N*T,1),0,M*N*T,M*N*T);                         R{5} = R{5}(r5,:);
L{5} = spdiags(ones((M-3)*N*T,1),0,(M-3)*N*T,M*N*T);                 L{5} = L{5}(:,l5);
K{5} = spdiags(ones((M-3)*(N-3)*T,1),0,(M-3)*(N-3)*T,(M-3)*(N-3)*T); K{5} = K{5}(:,k5);
% % I-K
% R{6} = spdiags(ones(M*N*T,1),0,M*N*T,M*N*T);                         R{6} = R{6}(r6,:);
% L{6} = spdiags(ones(M*(N-3)*T,1),0,M*(N-3)*T,M*N*T);                 L{6} = L{6}(:,l6);
% K{6} = spdiags(ones(M*(N-3)*(T-3),1),0,M*(N-3)*(T-3),M*(N-3)*(T-3)); K{6} = K{6}(:,k6);
% % J-K
% R{7} = spdiags(ones(M*N*T,1),0,M*N*T,M*N*T);                         R{7} = R{7}(r7,:);
% L{7} = spdiags(ones(M*N*(T-3),1),0,M*N*(T-3),M*N*(T-3));             L{7} = L{7}(:,l7);
% K{7} = spdiags(ones((M-3)*N*(T-3),1),0,(M-3)*N*(T-3),(M-3)*N*(T-3)); K{7} = K{7}(:,k7);

%% FOR ORDER 1 
% I-J
AA{1}{1} = K{5}.'*L{5}.'*R{5}*B{2}; % average along J for D1
AA{1}{2} = K{5}.'*L{5}.'*R{5}*B{1}; % average along I for D2
% I-K
AA{1}{3} = K{5}.'*L{5}.'*R{5}*B{3}; % average along K for D1
AA{1}{4} = K{5}.'*L{5}.'*R{5}*B{1}; % average along I for D3
% J-K
AA{1}{5} = K{5}.'*L{5}.'*R{5}*B{3}; % average along K for D2
AA{1}{6} = K{5}.'*L{5}.'*R{5}*B{2}; % average along J for D3

%% FOR ORDER 2
% I-J
AA{2}{1}  = K{5}.'*L{5}.'*R{5}*B{4};  % average along J for D1*D1
AA{2}{2}  = K{5}.'*L{5}.'*R{5};       % already on cell centres
AA{2}{3}  = K{5}.'*L{5}.'*R{5};       % already on cell centres
AA{2}{4}  = K{5}.'*L{5}.'*R{5}*B{5};  % average along I for D2*D2
% I-K
AA{2}{5}  = K{5}.'*L{5}.'*R{5}*B{6};  % average along K for D1*D1
AA{2}{6}  = K{5}.'*L{5}.'*R{5};       % already on cell centres
AA{2}{7}  = K{5}.'*L{5}.'*R{5};       % already on cell centres
AA{2}{8}  = K{5}.'*L{5}.'*R{5}*B{7};  % average along I for D3*D3
% J-K
AA{2}{9}  = K{5}.'*L{5}.'*R{5}*B{8};  % average along K for D2*D2
AA{2}{10}  = K{5}.'*L{5}.'*R{5};       % already on cell centres
AA{2}{11}  = K{5}.'*L{5}.'*R{5};       % already on cell centres
AA{2}{12} = K{5}.'*L{5}.'*R{5}*B{9}; % average along J for D3*D3

% %% 2 DERIVATIVES
% AA{2}{1} = L{1}.'*R{1}*B{3}; % average for D11
% AA{2}{2} = L{4}.'*R{4};      % THEY ARE ALREADY ON THE CELL CENTRE
% AA{2}{3} = L{4}.'*R{4};      % THEY ARE ALREADY ON THE CELL CENTRE
% AA{2}{4} = L{2}.'*R{2}*B{4}; % average for D22
% 
% %% 3 DERIVATIVES
% % this works
% AA{3}{1} = L{1}.'*R{1}*B{2};
% AA{3}{2} = L{1}.'*R{1}*B{1};
% AA{3}{3} = L{2}.'*R{2}*B{2};
% AA{3}{4} = L{2}.'*R{2}*B{2};
% AA{3}{5} = L{1}.'*R{1}*B{1};
% AA{3}{6} = L{2}.'*R{2}*B{1};
% AA{3}{1} = L{1}.'*R{1}*B{3};
% AA{3}{2} = L{1}.'*R{1}*B{1};
% AA{3}{3} = L{2}.'*R{2}*B{3};
% AA{3}{4} = L{2}.'*R{2}*B{3};
% AA{3}{5} = L{1}.'*R{1}*B{1};
% AA{3}{6} = L{2}.'*R{2}*B{1};



% % averaging for v
% C{1} = L{3}.'*R{3}*B{1};
% C{2} = L{3}.'*R{3}*B{2};