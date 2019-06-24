function Dbv_qa = discretization_order(u,A,params,D1,D2)

[M,N] = size(u);

D = discretization(u,params,D1,D2);

% A is the interpolation operator matrix that brings the gradient of u 
% in the correct position (i.e. cell centres)

switch params.order
    
    case 1     
        Dbv_qa(:,:,1) = reshape(A{1}{1}*flatten(D(:,:,1),1),M-3,N-3); % 1
        Dbv_qa(:,:,2) = reshape(A{1}{2}*flatten(D(:,:,2),1),M-3,N-3); % 2
    
    case 2
        Dbv_qa(:,:,1) = reshape(A{2}{1}*flatten(D(:,:,1),1),M-3,N-3); % 11
        Dbv_qa(:,:,2) = reshape(A{2}{2}*flatten(D(:,:,2),1),M-3,N-3); % 12
        Dbv_qa(:,:,3) = reshape(A{2}{3}*flatten(D(:,:,3),1),M-3,N-3); % 21
        Dbv_qa(:,:,4) = reshape(A{2}{4}*flatten(D(:,:,4),1),M-3,N-3); % 22
    
    case 3  
        Dbv_qa(:,:,1) = reshape(A{3}{1}*flatten(D(:,:,1),1),M-3,N-3); % 111
        Dbv_qa(:,:,2) = reshape(A{3}{2}*flatten(D(:,:,2),1),M-3,N-3); % 121
        Dbv_qa(:,:,3) = reshape(A{3}{2}*flatten(D(:,:,3),1),M-3,N-3); % 211
        Dbv_qa(:,:,4) = reshape(A{3}{3}*flatten(D(:,:,4),1),M-3,N-3); % 221
        
        Dbv_qa(:,:,5) = reshape(A{3}{5}*flatten(D(:,:,5),1),M-3,N-3); % 112
        Dbv_qa(:,:,6) = reshape(A{3}{3}*flatten(D(:,:,6),1),M-3,N-3); % 122
        Dbv_qa(:,:,7) = reshape(A{3}{3}*flatten(D(:,:,7),1),M-3,N-3); % 212
        Dbv_qa(:,:,8) = reshape(A{3}{6}*flatten(D(:,:,8),1),M-3,N-3); % 222
        
end

end