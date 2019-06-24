function [u,timecpu,PSNR_TV3D] = ROF2Dt_orders(u,uorig,params)

%% MIRROR BOUNADRY
u  = padarray(u,[1 1,1],'replicate','both');

[M,N,T] = size(u);

%% --- START TV3D

%% COMPUTE THE DERIVATIVE OPERATOR OF ORDER Q FOR EACH Q, WEIGHTED BY M
% Compute derivatives along (i,j,k)
D  = Kmat(u);

%% COMPUTE SADDLE-POINT OPERATORS AND PROX
% F:= (lamdba/eta)*TDV  is the vectorial soft thresholding.
K         = @(u,params)       cat(4,reshape(D{1}*u(:),M,N,T),reshape(D{2}*u(:),M,N,T),reshape(D{3}*u(:),M,N,T));
KS        = @(y,q)            reshape( D{1}.'*reshape(y(:,:,:,1),M*N*T,1) + D{2}.'*reshape(y(:,:,:,2),M*N*T,1) +  D{3}.'*reshape(y(:,:,:,3),M*N*T,1),M,N,T );
ProxFS    = @(y,sigma,params) y./max(1,repmat(norms(y,2,4)./params.lambda,1,1,1,size(y,4)));

% G = 0.5*eta*|| u - u^\diamond||_2^2
ProxG     = @(x,u,tau)( x+tau*params.eta*u )/( 1+tau*params.eta );

%% PRIMAL DUAL via Chambolle-Pock
options.niter  = params.niter;
tstart  = cputime;
u       = perform_primal_dual(uorig, u, params, K, KS, ProxFS, ProxG, options);
timecpu = cputime-tstart;

u = u(2:end-1,2:end-1,2:end-1);
PSNR_TV3D = psnr(uorig,u);

end

function x = perform_primal_dual(xorig, x, params, K, KS, ProxFS, ProxG, options)
% perform_primal_dual - primal-dual algorithm
%
%    [x,R] = perform_admm(x, K,  KS, ProxFS, ProxG, options);
%
%   Solves
%       min_x F(K*x) + G(x)
%   where F and G are convex proper functions with an easy to compute proximal operator,
%   and where K is a linear operator
%
%   Uses the Preconditioned Alternating direction method of multiplier (ADMM) method described in
%       Antonin Chambolle, Thomas Pock,
%       A first-order primal-dual algorithm for convex problems with applications to imaging,
%       Preprint CMAP-685
%
%   INPUTS:
%   ProxFS(y,sigma) computes Prox_{sigma*F^*}(y)
%   ProxG(x,tau) computes Prox_{tau*G}(x)
%   K(y) is a linear operator.
%   KS(y) compute K^*(y) the dual linear operator.
%   options.sigma and options.tau are the parameters of the
%       method, they shoudl satisfy sigma*tau*norm(K)^2<1
%   options.theta=1 for the ADMM, but can be set in [0,1].
%   options.verb=0 suppress display of progression.
%   options.niter is the number of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS:
%   x is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%   Copyright (c) 2010 Gabriel Peyre

options.null = 0;
niter  = getoptions(options, 'niter', 100);
theta  = getoptions(options, 'theta', 1);

%%%% ADMM parameters %%%%
sigma = getoptions(options, 'sigma', -1);
tau   = getoptions(options, 'tau', -1);

% INITIALIZATION
psnr_local = zeros(niter,1);
PDGAP      = zeros(niter,1);

if sigma<0 || tau<0
    rr = randn(size(x));
end

Q = numel(params.lambda);
index_Q = find(params.lambda);

paramslocal = cell(Q,1);
xstar       = cell(Q,1);
y           = cell(Q,1);
Lq          = cell(Q,1);

for q=index_Q
    
    paramslocal{q}        = params;
    paramslocal{q}.order  = params.order(q);
    paramslocal{q}.lambda = params.lambda(q);
    
    xstar{q} = 0;
    y{q}     = 0;
    Lq{q}    = 0;
    
%     if sigma<0 || tau<0
%         Lq{q} = compute_operator_norm(@(x) KS(K(x,paramslocal{q}),q),rr);
%         y{q}  = K(x,paramslocal{q});
%     end

end

% OPERATOR NORM
%L = max([8.^find(params.lambda)]);
L = max([(24).^find(params.lambda)]);

if params.acceleration
    
    %  ADMM
    tau   = 1/sqrt(L);
    sigma = 1./(tau*L);
    %gamma = 0.5;
    gamma = 0.35*params.eta;
    %gamma = 0.35*eta;
    
else
    sigma = 10;
    tau   = 0.9/(sigma*L);
end

if params.verbose_text
    fprintf('\n  ITER  |   PSNR  |    GAP    |  diff dual  \n');
    fprintf('----------------------------------------------\n');
end

xhat   = x;
xnoise = x;


I_MAX = 1;
psnr_start = 10*log10(I_MAX^2/mean((xorig(:)-reshape(x(2:end-1,2:end-1,2:end-1),[],1)).^2));
if params.verbose_text
    fprintf('   %02d   |  %2.2f  |  %2.2e |  %2.2e\n',0,psnr_start,0,Inf);
end

diff_dual = Inf;
for iter = 1:niter
    
    xold = x;
    
    % DUAL PROBLEM
    for q=index_Q
        yold{q}  = y{q};
        y{q}     = ProxFS( y{q} + sigma*K(xhat,paramslocal{q}), sigma, paramslocal{q});
        xstar{q} = KS(y{q},q);
    end
    
    
    if iter>1
        diff_dual_local = 0;
        for q=index_Q
            diff_dual_local = diff_dual_local + rms(y{q}(:) - yold{q}(:));
        end
        diff_dual = diff_dual_local/numel(index_Q);
    end
    
    % PRIMAL PROBLEM
    %x = ProxG(  x-tau*sum(cat(4,xstar{index_Q}),4), xnoise, tau, eta);
    x = ProxG(  x-tau*sum(cat(4,xstar{index_Q}),4), xnoise, tau);
    
    % EXTRAPOLATION
    xhat = x + theta * (x-xold);
    
    % ACCELERATION
    if params.acceleration
        theta = 1./sqrt(1+2*gamma*tau);
        tau   = theta*tau;
        sigma = sigma/theta;
    end
    
    
    %psnr_local(iter) = psnr(xorig(:,:,:),x(2:end-1,2:end-1,2:end-1));
    I_MAX = 1;
    psnr_local(iter) = 10*log10(I_MAX^2/mean((xorig(:)-reshape(x(2:end-1,2:end-1,2:end-1),[],1)).^2));
    
    if params.verbose_text && ~mod(iter,5)
        fprintf('   %02d   |  %2.2f  |  %2.2e |  %2.2e \n',iter,psnr_local(iter),PDGAP(iter),diff_dual)
    end
    
    if diff_dual<params.tolerance
       break 
    end
    
end
end

function [L,e] = compute_operator_norm(A,n)
% compute_operator_norm - compute operator norm
%
%   [L,e] = compute_operator_norm(A,n);
%
%   Copyright (c) 2010 Gabriel Peyre

if length(n)==1
    u = randn(n,1); u = u/norm(u);
else
    u = n;
    u = u/norm(u);
end
e = [];
for i=1:30
    v = A(u);
    e(end+1) = sum(u(:).*v(:));
    u = v/norm(v(:));
end
L = e(end);
end

function v = getoptions(options, name, v, mendatory)
% getoptions - retrieve options parameter
%
%   v = getoptions(options, 'entry', v0);
% is equivalent to the code:
%   if isfield(options, 'entry')
%       v = options.entry;
%   else
%       v = v0;
%   end
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<4
    mendatory = 0;
end

if isfield(options, name)
    v = eval(['options.' name ';']);
elseif mendatory
    error(['You have to provide options.' name '.']);
end
end

function gap = compute_pdgap(x,v,y,xstar,paramslocal,B,index_Q, F, G, K)

gap = 0;

for q = index_Q
    gap = gap + F(K(x,v,paramslocal{q},B),paramslocal{q}) + sum(sum( sum(-y{q}.*K(x,v,paramslocal{q},B),3) ));
end

gap = gap + G(x + sum(cat(3,xstar{index_Q}),3)) ;

end

function DM = divMq(y,M,D1,D2,D3,B)

Q = numel(M);

switch Q
    case 1
        
        z1 = (B{1}{1}*D1).'*flatten(M{1}(:,:,:,1,1,1).*y(:,:,:,1) + M{1}(:,:,:,2,1,1).*y(:,:,:,2),1);
        z2 = (B{1}{2}*D2).'*flatten(M{1}(:,:,:,1,2,1).*y(:,:,:,1) + M{1}(:,:,:,2,2,1).*y(:,:,:,2),1);
        z3 = (B{1}{3}*D1).'*flatten(M{1}(:,:,:,1,1,2).*y(:,:,:,3) + M{1}(:,:,:,2,1,2).*y(:,:,:,4),1);
        z4 = (B{1}{4}*D3).'*flatten(M{1}(:,:,:,1,2,2).*y(:,:,:,3) + M{1}(:,:,:,2,2,2).*y(:,:,:,4),1);
        z5 = (B{1}{5}*D2).'*flatten(M{1}(:,:,:,1,1,3).*y(:,:,:,5) + M{1}(:,:,:,2,1,3).*y(:,:,:,6),1);
        z6 = (B{1}{6}*D3).'*flatten(M{1}(:,:,:,1,2,3).*y(:,:,:,5) + M{1}(:,:,:,2,2,3).*y(:,:,:,6),1);
        
        DM = reshape(z1 + z2 + z3 + z4 + z5 + z6, size(y,1)+3,size(y,2)+3,size(y,3)+3 );
        
    case 2       
        % {i,j}
        zz1   = (B{2}{1}*D1) .'* flatten(M{2}(:,:,:,1,1,1).*y(:,:,:,1)  + M{2}(:,:,:,2,1,1).*y(:,:,:,2),1);
        zz2   = (B{2}{2}*D2) .'* flatten(M{2}(:,:,:,1,2,1).*y(:,:,:,1)  + M{2}(:,:,:,2,2,1).*y(:,:,:,2),1);
        zz3   = (B{2}{3}*D1) .'* flatten(M{2}(:,:,:,1,1,1).*y(:,:,:,3)  + M{2}(:,:,:,2,1,1).*y(:,:,:,4),1);
        zz4   = (B{2}{4}*D2) .'* flatten(M{2}(:,:,:,1,2,1).*y(:,:,:,3)  + M{2}(:,:,:,2,2,1).*y(:,:,:,4),1);
         
        % {i,k}
        zz5   = (B{2}{5}*D1) .'* flatten(M{2}(:,:,:,1,1,2).*y(:,:,:,5)  + M{2}(:,:,:,2,1,2).*y(:,:,:,6),1);
        zz6   = (B{2}{6}*D3) .'* flatten(M{2}(:,:,:,1,2,2).*y(:,:,:,5)  + M{2}(:,:,:,2,2,2).*y(:,:,:,6),1);
        zz7   = (B{2}{7}*D1) .'* flatten(M{2}(:,:,:,1,1,2).*y(:,:,:,7)  + M{2}(:,:,:,2,1,2).*y(:,:,:,8),1);
        zz8   = (B{2}{8}*D3) .'* flatten(M{2}(:,:,:,1,2,2).*y(:,:,:,7)  + M{2}(:,:,:,2,2,2).*y(:,:,:,8),1);
       
        % {j,k}
        zz9   = (B{2}{9}*D2)  .'* flatten(M{2}(:,:,:,1,1,3).*y(:,:,:,9)  + M{2}(:,:,:,2,1,3).*y(:,:,:,10),1);
        zz10  = (B{2}{10}*D3) .'* flatten(M{2}(:,:,:,1,2,3).*y(:,:,:,9)  + M{2}(:,:,:,2,2,3).*y(:,:,:,10),1);
        zz11  = (B{2}{11}*D2) .'* flatten(M{2}(:,:,:,1,1,3).*y(:,:,:,11) + M{2}(:,:,:,2,1,3).*y(:,:,:,12),1);
        zz12  = (B{2}{12}*D3) .'* flatten(M{2}(:,:,:,1,2,3).*y(:,:,:,11) + M{2}(:,:,:,2,2,3).*y(:,:,:,12),1);
        
        % div
        ww1 = (zz1 + zz2);
        ww2 = (zz3 + zz4);        
        ww3 = (zz5 + zz6);
        ww4 = (zz7 + zz8);
        ww5 = (zz9 + zz10);
        ww6 = (zz11 + zz12);
        
        W1  = reshape(cat(2,ww1,ww2),size(M{2},1)+3,size(M{2},2)+3,size(M{2},3)+3,2);
        W2  = reshape(cat(2,ww3,ww4),size(M{2},1)+3,size(M{2},2)+3,size(M{2},3)+3,2);
        W3  = reshape(cat(2,ww5,ww6),size(M{2},1)+3,size(M{2},2)+3,size(M{2},3)+3,2);
        
        w1 = (B{1}{1}*D1) .'* flatten( M{1}(:,:,:,1,1,1).*W1(:,:,:,1)  + M{1}(:,:,:,2,1,1).*W1(:,:,:,2),1);
        w2 = (B{1}{2}*D2) .'* flatten( M{1}(:,:,:,1,2,1).*W1(:,:,:,1)  + M{1}(:,:,:,2,2,1).*W1(:,:,:,2),1);
        w3 = (B{1}{3}*D1) .'* flatten( M{1}(:,:,:,1,1,2).*W2(:,:,:,1)  + M{1}(:,:,:,2,1,2).*W2(:,:,:,2),1);
        w4 = (B{1}{4}*D3) .'* flatten( M{1}(:,:,:,1,2,2).*W2(:,:,:,1)  + M{1}(:,:,:,2,2,2).*W2(:,:,:,2),1);         
        w5 = (B{1}{5}*D2) .'* flatten( M{1}(:,:,:,1,1,3).*W3(:,:,:,1)  + M{1}(:,:,:,2,1,3).*W3(:,:,:,2),1);
        w6 = (B{1}{6}*D3) .'* flatten( M{1}(:,:,:,1,2,3).*W3(:,:,:,1)  + M{1}(:,:,:,2,2,3).*W3(:,:,:,2),1);
  
        % div
        z1 = sum(reshape(cat(2,w1,w2),size(M{1},1),size(M{1},2),size(M{1},3),2),4);
        z2 = sum(reshape(cat(2,w3,w4),size(M{1},1),size(M{1},2),size(M{1},3),2),4);
        z3 = sum(reshape(cat(2,w5,w6),size(M{1},1),size(M{1},2),size(M{1},3),2),4);   
        
        DM = z1+z2+z3;
end
    
end