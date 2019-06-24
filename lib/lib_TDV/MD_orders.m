function Z = MD_orders(u,LR,B,D1,D2,D3)

Q       = numel(LR);
[M,N,T] = size(u);

switch Q
    case 0
        % u
        Z = u;
        
    case 1
        % on each color channe do LR*grad
        %Z = basicMD(LR{1},B,D1,D2,D3);
        [Z{1},Z{2}] = basicMD(LR{1}(:,:,:,:,:,1),B{1}{1}*D1,B{1}{2}*D2);
        [Z{3},Z{4}] = basicMD(LR{1}(:,:,:,:,:,2),B{1}{3}*D1,B{1}{4}*D3);
        [Z{5},Z{6}] = basicMD(LR{1}(:,:,:,:,:,3),B{1}{5}*D2,B{1}{6}*D3);
        
    case 2
     
        % MD
        %MD = basicMD(LR{1},{B{1}},D1,D2,D3);
        [MD{1},MD{2}] = basicMD(LR{1}(:,:,:,:,:,1),B{1}{1}*D1,B{1}{2}*D2);
        [MD{3},MD{4}] = basicMD(LR{1}(:,:,:,:,:,2),B{1}{3}*D1,B{1}{4}*D3);
        [MD{5},MD{6}] = basicMD(LR{1}(:,:,:,:,:,3),B{1}{5}*D2,B{1}{6}*D3);
        
        % MD_MD
        % {i,j}
        [Z{1},Z{2}]   = basicMD(LR{2}(:,:,:,:,:,1), B{2}{1}*D1*MD{1}, B{2}{2}*D2*MD{1});
        [Z{3},Z{4}]   = basicMD(LR{2}(:,:,:,:,:,1), B{2}{3}*D1*MD{2}, B{2}{4}*D2*MD{2});
        % {i,k}
        [Z{5},Z{6}]   = basicMD(LR{2}(:,:,:,:,:,2), B{2}{5}*D1*MD{3}, B{2}{6}*D3*MD{3});
        [Z{7},Z{8}]   = basicMD(LR{2}(:,:,:,:,:,2), B{2}{7}*D1*MD{4}, B{2}{8}*D3*MD{4});
        % {j,k}
        [Z{9},Z{10}]  = basicMD(LR{2}(:,:,:,:,:,3), B{2}{9}*D2*MD{5}, B{2}{10}*D3*MD{5});
        [Z{11},Z{12}] = basicMD(LR{2}(:,:,:,:,:,3), B{2}{11}*D2*MD{6},B{2}{12}*D3*MD{6});
        
    case 3   
        
        % MD
        [MD{1},MD{2}] = basicMD(LR{1}(:,:,:,:,:,1),D1,D2);
        [MD{3},MD{4}] = basicMD(LR{1}(:,:,:,:,:,2),D1,D3);
        [MD{5},MD{6}] = basicMD(LR{1}(:,:,:,:,:,3),D2,D3);
        
        % MD_MD
        % {i,j}
        [MD_MD{1},MD_MD{2}]   = basicMD(LR{2}(:,:,:,:,:,1), B{2}{1}*D1*MD{1}, B{2}{2}*D2*MD{1});
        [MD_MD{3},MD_MD{4}]   = basicMD(LR{2}(:,:,:,:,:,1), B{2}{3}*D1*MD{2}, B{2}{4}*D2*MD{2});
        % {i,k}
        [MD_MD{5},MD_MD{6}]   = basicMD(LR{2}(:,:,:,:,:,2), B{2}{5}*D1*MD{3}, B{2}{6}*D3*MD{3});
        [MD_MD{7},MD_MD{8}]   = basicMD(LR{2}(:,:,:,:,:,2), B{2}{7}*D1*MD{4}, B{2}{8}*D3*MD{4});
        % {j,k}
        [MD_MD{9},MD_MD{10}]  = basicMD(LR{2}(:,:,:,:,:,3), B{2}{9}*D2*MD{5}, B{2}{10}*D3*MD{5});
        [MD_MD{11},MD_MD{12}] = basicMD(LR{2}(:,:,:,:,:,3), B{2}{11}*D2*MD{6},B{2}{12}*D3*MD{6});
        
        % MD_MD_MD
        % {i,j}
        [Z{1}{1}{1},Z{2}{1}{1}] = basicMD(LR{3}(:,:,:,:,:,1),B{3}{1}*D1*MD_MD{1}{1},B{3}{5}*D2*MD_MD{1}{1});
        [Z{1}{1}{2},Z{2}{1}{2}] = basicMD(LR{3}(:,:,:,:,:,1),B{3}{2}*D1*MD_MD{1}{2},B{3}{3}*D2*MD_MD{1}{2});
        [Z{1}{2}{1},Z{2}{2}{1}] = basicMD(LR{3}(:,:,:,:,:,1),B{3}{2}*D1*MD_MD{2}{1},B{3}{3}*D2*MD_MD{2}{1});
        [Z{1}{2}{2},Z{2}{2}{2}] = basicMD(LR{3}(:,:,:,:,:,1),B{3}{3}*D1*MD_MD{2}{2},B{3}{6}*D2*MD_MD{2}{2});
        
        % {i,k} - B{3}{1} ETC TO BE CHANGED
        [Z{3}{3}{3},Z{4}{3}{3}] = basicMD(LR{3}(:,:,:,:,:,2),B{3}{1}*D1*MD_MD{5}{1},B{3}{5}*D3*MD_MD{1}{1});
        [Z{3}{3}{4},Z{4}{3}{4}] = basicMD(LR{3}(:,:,:,:,:,2),B{3}{2}*D1*MD_MD{1}{2},B{3}{3}*D3*MD_MD{1}{2});
        [Z{3}{4}{3},Z{4}{4}{3}] = basicMD(LR{3}(:,:,:,:,:,2),B{3}{2}*D1*MD_MD{2}{1},B{3}{3}*D3*MD_MD{2}{1});
        [Z{3}{4}{4},Z{4}{4}{4}] = basicMD(LR{3}(:,:,:,:,:,2),B{3}{3}*D1*MD_MD{2}{2},B{3}{6}*D3*MD_MD{2}{2});
        
        % {j,k} - B{3}{1} ETC TO BE CHANGED
        [Z{5}{5}{5},Z{6}{5}{5}] = basicMD(LR{3}(:,:,:,:,:,3),B{3}{1}*D2*MD_MD{1}{1},B{3}{5}*D3*MD_MD{1}{1});
        [Z{5}{5}{6},Z{6}{5}{6}] = basicMD(LR{3}(:,:,:,:,:,3),B{3}{2}*D2*MD_MD{1}{2},B{3}{3}*D3*MD_MD{1}{2});
        [Z{5}{6}{5},Z{6}{6}{5}] = basicMD(LR{3}(:,:,:,:,:,3),B{3}{2}*D2*MD_MD{2}{1},B{3}{3}*D3*MD_MD{2}{1});
        [Z{5}{6}{6},Z{6}{6}{6}] = basicMD(LR{3}(:,:,:,:,:,3),B{3}{3}*D2*MD_MD{2}{2},B{3}{6}*D3*MD_MD{2}{2});
end

return

% function Z = basicMD(LR,B,D1,D2,D3)
% 
% M = size(LR,1);
% N = size(LR,2);
% T = size(LR,3);
% 
% % {i,j}
% Z{1} = spdiags(flatten(LR(:,:,:,1,1,1),1),0,M*N*T,M*N*T)*B{1}{1}*D1 +...
%     spdiags(flatten(LR(:,:,:,1,2,1),1),0,M*N*T,M*N*T)*B{1}{2}*D2;
% 
% Z{2} = spdiags(flatten(LR(:,:,:,2,1,1),1),0,M*N*T,M*N*T)*B{1}{1}*D1 +...
%     spdiags(flatten(LR(:,:,:,2,2,1),1),0,M*N*T,M*N*T)*B{1}{2}*D2;
% 
% 
% % {i,k}
% Z{3} = spdiags(flatten(LR(:,:,:,1,1,2),1),0,M*N*T,M*N*T)*B{1}{3}*D1 +...
%     spdiags(flatten(LR(:,:,:,1,2,2),1),0,M*N*T,M*N*T)*B{1}{4}*D3;
% 
% Z{4} = spdiags(flatten(LR(:,:,:,2,1,2),1),0,M*N*T,M*N*T)*B{1}{3}*D1 +...
%     spdiags(flatten(LR(:,:,:,2,2,2),1),0,M*N*T,M*N*T)*B{1}{4}*D3;
% 
% % {j,k}
% Z{5} = spdiags(flatten(LR(:,:,:,1,1,3),1),0,M*N*T,M*N*T)*B{1}{5}*D2 +...
%     spdiags(flatten(LR(:,:,:,1,2,3),1),0,M*N*T,M*N*T)*B{1}{6}*D3;
% 
% Z{6} = spdiags(flatten(LR(:,:,:,2,1,3),1),0,M*N*T,M*N*T)*B{1}{5}*D2 +...
%     spdiags(flatten(LR(:,:,:,2,2,3),1),0,M*N*T,M*N*T)*B{1}{6}*D3;
% 
% 
% return

function [Z1,Z2] = basicMD(LR,D1,D2)

M = size(LR,1);
N = size(LR,2);
T = size(LR,3);

% {i,j}
Z1 = spdiags(flatten(LR(:,:,:,1,1),1),0,M*N*T,M*N*T)*D1 +...
    spdiags(flatten(LR(:,:,:,1,2),1),0,M*N*T,M*N*T)*D2;

Z2 = spdiags(flatten(LR(:,:,:,2,1),1),0,M*N*T,M*N*T)*D1 +...
    spdiags(flatten(LR(:,:,:,2,2),1),0,M*N*T,M*N*T)*D2;


return