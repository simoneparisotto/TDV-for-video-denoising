function [V,LAMBDA,A] = compute_structure_Tensor_sliced(u,sigma,rho,location)

[M,N,T] = size(u(2:end-1,2:end-1,2:end-1));

%% MEMORY STORAGE
lambda = zeros(M+1,N+1,T+1,3);   % on extended cell centres
v      = zeros(M+1,N+1,T+1,6);   % on extended cell centres
LAMBDA = zeros(M-1,N-1,T-1,6);   % on interior cell centres
V      = zeros(M-1,N-1,T-1,6);   % on interior cell centres

%% COMPUTE DERIVATIVES
D     = Kmat(u);

%% SMOOTHING PARAMETERS
epsilon       = 1e-15;
%filter_domain = 'frequency';
filter_domain = 'spatial';
fs   = 15;

%% PRE-SMOOTHING: noise scale
usigma = imgaussfilt3(u,sigma,'FilterDomain',filter_domain,'FilterSize',fs);
%usigma = smooth3(u,'gaussian',fs,sigma);

%% COMPUTE DERIVATIVES ALONG I-J-K on v_extra
DU = cat(4,interpn(location.ii_u1_extra,location.jj_u1_extra,location.kk_u1_extra, reshape(D{1}*usigma(:),M+2,N+2,T+2), location.ii_v_extra,location.jj_v_extra,location.kk_v_extra),...
           interpn(location.ii_u2_extra,location.jj_u2_extra,location.kk_u2_extra, reshape(D{2}*usigma(:),M+2,N+2,T+2), location.ii_v_extra,location.jj_v_extra,location.kk_v_extra),...
           interpn(location.ii_u3_extra,location.jj_u3_extra,location.kk_u3_extra, reshape(D{3}*usigma(:),M+2,N+2,T+2), location.ii_v_extra,location.jj_v_extra,location.kk_v_extra));

%% COMPUTE TENSOR AND POST-SMOOTHING VIA RHO
DUrho = zeros(M+1,N+1,T+1,3,3);
DUrho(:,:,:,1,1) = imgaussfilt3(DU(:,:,:,1).*DU(:,:,:,1),rho,'FilterDomain',filter_domain,'FilterSize',fs);
DUrho(:,:,:,1,2) = imgaussfilt3(DU(:,:,:,1).*DU(:,:,:,2),rho,'FilterDomain',filter_domain,'FilterSize',fs);
DUrho(:,:,:,1,3) = imgaussfilt3(DU(:,:,:,1).*DU(:,:,:,3),rho,'FilterDomain',filter_domain,'FilterSize',fs);
DUrho(:,:,:,2,1) = imgaussfilt3(DU(:,:,:,2).*DU(:,:,:,1),rho,'FilterDomain',filter_domain,'FilterSize',fs);
DUrho(:,:,:,2,2) = imgaussfilt3(DU(:,:,:,2).*DU(:,:,:,2),rho,'FilterDomain',filter_domain,'FilterSize',fs);
DUrho(:,:,:,2,3) = imgaussfilt3(DU(:,:,:,2).*DU(:,:,:,3),rho,'FilterDomain',filter_domain,'FilterSize',fs);
DUrho(:,:,:,3,1) = imgaussfilt3(DU(:,:,:,3).*DU(:,:,:,1),rho,'FilterDomain',filter_domain,'FilterSize',fs);
DUrho(:,:,:,3,2) = imgaussfilt3(DU(:,:,:,3).*DU(:,:,:,2),rho,'FilterDomain',filter_domain,'FilterSize',fs);
DUrho(:,:,:,3,3) = imgaussfilt3(DU(:,:,:,3).*DU(:,:,:,3),rho,'FilterDomain',filter_domain,'FilterSize',fs);
% DUrho(:,:,:,1,1) = smooth3(DU(:,:,:,1).*DU(:,:,:,1),'gaussian',fs,rho);
% DUrho(:,:,:,1,2) = smooth3(DU(:,:,:,1).*DU(:,:,:,2),'gaussian',fs,rho);
% DUrho(:,:,:,1,3) = smooth3(DU(:,:,:,1).*DU(:,:,:,3),'gaussian',fs,rho);
% DUrho(:,:,:,2,1) = smooth3(DU(:,:,:,2).*DU(:,:,:,1),'gaussian',fs,rho);
% DUrho(:,:,:,2,2) = smooth3(DU(:,:,:,2).*DU(:,:,:,2),'gaussian',fs,rho);
% DUrho(:,:,:,2,3) = smooth3(DU(:,:,:,2).*DU(:,:,:,3),'gaussian',fs,rho);
% DUrho(:,:,:,3,1) = smooth3(DU(:,:,:,3).*DU(:,:,:,1),'gaussian',fs,rho);
% DUrho(:,:,:,3,2) = smooth3(DU(:,:,:,3).*DU(:,:,:,2),'gaussian',fs,rho);
% DUrho(:,:,:,3,3) = smooth3(DU(:,:,:,3).*DU(:,:,:,3),'gaussian',fs,rho);


delta     = @(S)(S(:,:,1)-S(:,:,2)).^2 + 4*S(:,:,3).^2;
eigenval  = @(S)deal( ...
    (S(:,:,1)+S(:,:,2)+sqrt(delta(S)))/2,  ...
    (S(:,:,1)+S(:,:,2)-sqrt(delta(S)))/2 );
normalize = @(u)u./repmat(sqrt(sum(u.^2,3)), [1 1 2]);
eig1      = @(S)normalize( cat(3,2*S(:,:,3), S(:,:,2)-S(:,:,1)+sqrt(delta(S)) ) );
ortho     = @(u) cat(3,-u(:,:,2), u(:,:,1));
eigbasis  = @(S) (eig1(S));

% I-J
DUrho_IJ = cat(4,DUrho(:,:,:,1,1),DUrho(:,:,:,2,2),DUrho(:,:,:,1,2));
for kk=1:T+1
    [lambda(:,:,kk,1),lambda(:,:,kk,2)] = eigenval(squeeze(DUrho_IJ(:,:,kk,:)));
    v(:,:,kk,[1,2])                     = eigbasis(squeeze(DUrho_IJ(:,:,kk,:)));
end
% I-K
DUrho_IK = cat(4,DUrho(:,:,:,1,1),DUrho(:,:,:,3,3),DUrho(:,:,:,1,3));
for jj=1:N+1
    [lambda(:,jj,:,3),lambda(:,jj,:,4)] = eigenval(squeeze(DUrho_IK(:,jj,:,:)));
    v(:,jj,:,[3,4])                     = eigbasis(squeeze(DUrho_IK(:,jj,:,:)));
end

% J-K
DUrho_JK = cat(4,DUrho(:,:,:,2,2),DUrho(:,:,:,3,3),DUrho(:,:,:,2,3));
for ii=1:M+1
    [lambda(ii,:,:,5),lambda(ii,:,:,6)] = eigenval(squeeze(DUrho_JK(ii,:,:,:)));
    v(ii,:,:,[5,6])                     = eigbasis(squeeze(DUrho_JK(ii,:,:,:)));
end

% lambda(:,:,:,2) = lambda(:,:,:,1) + (1-lambda(:,:,:,1)).*exp(-1./(lambda(:,:,:,1)-lambda(:,:,:,2)).^2);
% lambda(:,:,:,4) = lambda(:,:,:,3) + (1-lambda(:,:,:,3)).*exp(-1./(lambda(:,:,:,3)-lambda(:,:,:,4)).^2);
% lambda(:,:,:,6) = lambda(:,:,:,5) + (1-lambda(:,:,:,5)).*exp(-1./(lambda(:,:,:,5)-lambda(:,:,:,6)).^2);

%% COMPUTE EIGENVALUES LAMBDA AND EIGENVECTORS V on the interior of the grid
% {i,j}
LAMBDA(:,:,:,1) = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  lambda(:,:,:,1), location.ii_v,location.jj_v,location.kk_v);
LAMBDA(:,:,:,2) = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  lambda(:,:,:,2), location.ii_v,location.jj_v,location.kk_v);
V(:,:,:,1)      = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  v(:,:,:,1), location.ii_v,location.jj_v,location.kk_v);
V(:,:,:,2)      = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  v(:,:,:,2), location.ii_v,location.jj_v,location.kk_v);
V(:,:,:,1:2)    = V(:,:,:,1:2) ./ repmat(sum(V(:,:,:,1:2).^2+epsilon,4).^0.5,[1,1,1,2]);
% {i,k}
LAMBDA(:,:,:,3) = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  lambda(:,:,:,3), location.ii_v,location.jj_v,location.kk_v);
LAMBDA(:,:,:,4) = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  lambda(:,:,:,4), location.ii_v,location.jj_v,location.kk_v);
V(:,:,:,3)      = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  v(:,:,:,3), location.ii_v,location.jj_v,location.kk_v);
V(:,:,:,4)      = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  v(:,:,:,4), location.ii_v,location.jj_v,location.kk_v);
V(:,:,:,3:4)    = V(:,:,:,3:4) ./ repmat(sum(V(:,:,:,3:4).^2+epsilon,4).^0.5,[1,1,1,2]);
% {j,k}
LAMBDA(:,:,:,5) = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  lambda(:,:,:,5), location.ii_v,location.jj_v,location.kk_v);
LAMBDA(:,:,:,6) = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  lambda(:,:,:,6), location.ii_v,location.jj_v,location.kk_v);
V(:,:,:,5)      = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  v(:,:,:,5), location.ii_v,location.jj_v,location.kk_v);
V(:,:,:,6)      = interpn(location.ii_v_extra,location.jj_v_extra,location.kk_v_extra,  v(:,:,:,6), location.ii_v,location.jj_v,location.kk_v);
V(:,:,:,5:6)    = V(:,:,:,5:6) ./ repmat(sum(V(:,:,:,5:6).^2+epsilon,4).^0.5,[1,1,1,2]);

%% ANISOTROPY
%A(:,:,:,1)  = 1 - (abs(LAMBDA(:,:,:,1) - LAMBDA(:,:,:,2)) ./ ( epsilon + LAMBDA(:,:,:,1) + LAMBDA(:,:,:,2)) );
%A(:,:,:,2)  = 1 - (abs(LAMBDA(:,:,:,3) - LAMBDA(:,:,:,4)) ./ ( epsilon + LAMBDA(:,:,:,3) + LAMBDA(:,:,:,4)) );
%A(:,:,:,3)  = 1 - (abs(LAMBDA(:,:,:,5) - LAMBDA(:,:,:,6)) ./ ( epsilon + LAMBDA(:,:,:,5) + LAMBDA(:,:,:,6)) );

%A(:,:,:,1)  = 1-(abs(LAMBDA(:,:,:,1) - LAMBDA(:,:,:,2))./(epsilon + LAMBDA(:,:,:,1) + LAMBDA(:,:,:,2))).^2;
%A(:,:,:,2)  = 1-(abs(LAMBDA(:,:,:,3) - LAMBDA(:,:,:,4))./(epsilon + LAMBDA(:,:,:,3) + LAMBDA(:,:,:,4))).^2;
%A(:,:,:,3)  = 1-(abs(LAMBDA(:,:,:,5) - LAMBDA(:,:,:,6))./(epsilon + LAMBDA(:,:,:,5) + LAMBDA(:,:,:,6))).^2;


% A(:,:,:,1)   = abs(LAMBDA(:,:,:,2)./(epsilon + LAMBDA(:,:,:,1)));
% A(:,:,:,2)   = abs(LAMBDA(:,:,:,4)./(epsilon + LAMBDA(:,:,:,3)));
% A(:,:,:,3)   = abs(LAMBDA(:,:,:,6)./(epsilon + LAMBDA(:,:,:,5)));

A(:,:,:,1)   = (abs(LAMBDA(:,:,:,2)./(epsilon + LAMBDA(:,:,:,1))));
A(:,:,:,2)   = ones(size(LAMBDA(:,:,:,2)));
A(:,:,:,3)   = (abs(LAMBDA(:,:,:,4)./(epsilon + LAMBDA(:,:,:,3))));
A(:,:,:,4)   = ones(size(LAMBDA(:,:,:,3)));
A(:,:,:,5)   = (abs(LAMBDA(:,:,:,6)./(epsilon + LAMBDA(:,:,:,5))));
A(:,:,:,6)   = ones(size(LAMBDA(:,:,:,5)));


A(:,:,:,1)   = sqrt(abs(LAMBDA(:,:,:,2)./(epsilon + LAMBDA(:,:,:,1))));
A(:,:,:,2)   = ones(size(LAMBDA(:,:,:,2)));
A(:,:,:,3)   = sqrt(abs(LAMBDA(:,:,:,4)./(epsilon + LAMBDA(:,:,:,3))));
A(:,:,:,4)   = ones(size(LAMBDA(:,:,:,3)));
A(:,:,:,5)   = sqrt(abs(LAMBDA(:,:,:,6)./(epsilon + LAMBDA(:,:,:,5))));
A(:,:,:,6)   = ones(size(LAMBDA(:,:,:,5)));


