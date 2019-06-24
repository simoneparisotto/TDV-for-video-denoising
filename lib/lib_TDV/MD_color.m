function Z = MD_color(u,LR,B,D1,D2,D3,location)

Q         = numel(LR);
[M,N,C,T] = size(u);

switch Q
    case 0
        % u
        Z = u;
        
    case 1
        % on each color channe do LR*grad
        Z = basicMD(LR{1},B,D1,D2,D3,M,N,C,T);
        
    case 2
     
        % MD
        [MD{1},MD{2}] = basicMD(M{1},B{1}{1}*D1,B{1}{2}*D2,size(M{1},1),size(M{1},2));
        
        % MD_MD
        [Z{1}{1},Z{2}{1}] =  basicMD(M{2},B{2}{1}*D1*MD{1},B{2}{2}*D2*MD{1},size(M{2},1),size(M{2},2));
        [Z{1}{2},Z{2}{2}] =  basicMD(M{2},B{2}{3}*D1*MD{2},B{2}{4}*D2*MD{2},size(M{2},1),size(M{2},2));
        
    case 3   
        
        % MD
        [MD{1},MD{2}]    = basicMD(M{1},B{1}{1}*D1,B{1}{1}*D2,size(M{1},1),size(M{1},2));
        
        % MD_MD
        [MD_MD{1}{1},MD_MD{2}{1}] = basicMD(M{2},B{2}{1}*D1*MD{1},B{2}{2}*D2*MD{1},size(M{2},1),size(M{2},2));
        [MD_MD{1}{2},MD_MD{2}{2}] = basicMD(M{2},B{2}{3}*D1*MD{2},B{2}{4}*D2*MD{2},size(M{2},1),size(M{2},2));

        % MD_MD_MD
        [Z{1}{1}{1},Z{2}{1}{1}] = basicMD(M{3},B{3}{1}*D1*MD_MD{1}{1},B{3}{5}*D2*MD_MD{1}{1},size(M{3},1),size(M{3},2));
        [Z{1}{1}{2},Z{2}{1}{2}] = basicMD(M{3},B{3}{2}*D1*MD_MD{1}{2},B{3}{3}*D2*MD_MD{1}{2},size(M{3},1),size(M{3},2));
        [Z{1}{2}{1},Z{2}{2}{1}] = basicMD(M{3},B{3}{2}*D1*MD_MD{2}{1},B{3}{3}*D2*MD_MD{2}{1},size(M{3},1),size(M{3},2));
        [Z{1}{2}{2},Z{2}{2}{2}] = basicMD(M{3},B{3}{3}*D1*MD_MD{2}{2},B{3}{6}*D2*MD_MD{2}{2},size(M{3},1),size(M{3},2));
end

return

function Z = basicMD(LR,B,D1,D2,D3,M,N,C,T)

M = M-3;
N = N-3;
T = T-3;

for cc=1:C
    
    LRloc = squeeze(LR(:,:,cc,:,:,:,:));
    
    % {i,j}
    Z{cc}{1} = spdiags(flatten(LRloc(:,:,:,1,1,1),1),0,M*N*T,M*N*T)*B{1}{1}*D1 +...
               spdiags(flatten(LRloc(:,:,:,1,2,1),1),0,M*N*T,M*N*T)*B{1}{2}*D2;
    
    Z{cc}{2} = spdiags(flatten(LRloc(:,:,:,2,1,1),1),0,M*N*T,M*N*T)*B{1}{1}*D1 +...
               spdiags(flatten(LRloc(:,:,:,2,2,1),1),0,M*N*T,M*N*T)*B{1}{2}*D2;
    
    
    % {i,k}
    Z{cc}{3} = spdiags(flatten(LRloc(:,:,:,1,1,2),1),0,M*N*T,M*N*T)*B{1}{3}*D1 +...
               spdiags(flatten(LRloc(:,:,:,1,2,2),1),0,M*N*T,M*N*T)*B{1}{4}*D3;
    
    Z{cc}{4} = spdiags(flatten(LRloc(:,:,:,2,1,2),1),0,M*N*T,M*N*T)*B{1}{3}*D1 +...
               spdiags(flatten(LRloc(:,:,:,2,2,2),1),0,M*N*T,M*N*T)*B{1}{4}*D3;
    
    % {j,k}
    Z{cc}{5} = spdiags(flatten(LRloc(:,:,:,1,1,3),1),0,M*N*T,M*N*T)*B{1}{5}*D2 +...
               spdiags(flatten(LRloc(:,:,:,1,2,3),1),0,M*N*T,M*N*T)*B{1}{6}*D3;
    
    Z{cc}{6} = spdiags(flatten(LRloc(:,:,:,2,1,3),1),0,M*N*T,M*N*T)*B{1}{5}*D2 +...
               spdiags(flatten(LRloc(:,:,:,2,2,3),1),0,M*N*T,M*N*T)*B{1}{6}*D3;
    
end

return