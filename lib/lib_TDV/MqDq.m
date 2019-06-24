function w = MqDq(u,M,B,D1,D2,location)

Q = numel(M);

switch Q
    case 0
        % u
        w = u;
    case 1
        % MDu
        [w1,w2] = MD(u,M{1},B,D1,D2,location);
        w = cat(3,w1,w2);
    case 2
        % MDu
        [w1,w2]   = MD(u,M{1},B,D1,D2,location);
        % MD(MDu)
        [w11,w12] = MD(w1,M{2},B,D1,D2,location);
        [w21,w22] = MD(w2,M{2},B,D1,D2,location);
        w = cat(3,w11,w12,w21,w22);
    case 3
        % MDu
        [w1,w2]  = MD(u,M{1},B,D1,D2,location);
        % MD(MDu)
        [w11,w12] = MD(w1,M{2},B,D1,D2,location);
        [w21,w22] = MD(w2,M{2},B,D1,D2,location);
        % MD(MD(MDu))
        [w111,w112] = MD(w11,M{3},B,D1,D2,location);
        [w121,w122] = MD(w12,M{3},B,D1,D2,location);
        [w211,w212] = MD(w21,M{3},B,D1,D2,location);
        [w221,w222] = MD(w22,M{3},B,D1,D2,location);
        w = cat(3,w111,w112,w121,w122,w211,w212,w221,w222);
end

return

function [MDU1,MDU2] = MD(y,M,B,D1,D2,location)
    DU{1} = reshape(D1*flatten(y,1),size(M,1),size(M,2));
    DU{2} = reshape(D2*flatten(y,1),size(M,1),size(M,2));
    
    MDU1 = M(:,:,1).*DU{1} + M(:,:,2).*DU{2};
    MDU2 = M(:,:,3).*DU{1} + M(:,:,4).*DU{2};

return


% function w = MqDq(u,M,B,D1,D2,location)
% 
% Q = numel(M);
% 
% switch Q
%     case 0
%         % u
%         w = u;
%     case 1
%         % MDu
%         [w1,w2] = MD(u,M{1},B,D1,D2,location);
%         w = cat(3,w1,w2);
%     case 2
%         % MDu
%         [w1,w2]   = MD(u,M{1},B,D1,D2,location);
%         % MD(MDu)
%         [w11,w12] = MD(w1,M{2},B,D1,D2,location);
%         [w21,w22] = MD(w2,M{2},B,D1,D2,location);
%         w = cat(3,w11,w12,w21,w22);
%     case 3
%         % MDu
%         [w1,w2]  = MD(u,M{1},B,D1,D2,location);
%         % MD(MDu)
%         w1 = interpn(location.ii_v,location.jj_v,w1,location.ii_u_extra,location.jj_u_extra,'linear',0);
%         w2 = interpn(location.ii_v,location.jj_v,w2,location.ii_u_extra,location.jj_u_extra,'linear',0);
%         [w11,w12] = MD(w1,M{2},B,D1,D2,location);
%         [w21,w22] = MD(w2,M{2},B,D1,D2,location);
%         % MD(MD(MDu))
%         w11 = interpn(location.ii_v,location.jj_v,w11,location.ii_u_extra,location.jj_u_extra,'linear',0);
%         w12 = interpn(location.ii_v,location.jj_v,w12,location.ii_u_extra,location.jj_u_extra,'linear',0);
%         w21 = interpn(location.ii_v,location.jj_v,w21,location.ii_u_extra,location.jj_u_extra,'linear',0);
%         w22 = interpn(location.ii_v,location.jj_v,w22,location.ii_u_extra,location.jj_u_extra,'linear',0);
%         [w111,w112] = MD(w11,M{3},B,D1,D2,location);
%         [w121,w122] = MD(w12,M{3},B,D1,D2,location);
%         [w211,w212] = MD(w21,M{3},B,D1,D2,location);
%         [w221,w222] = MD(w22,M{3},B,D1,D2,location);
%         w = cat(3,w111,w112,w121,w122,w211,w212,w221,w222);
% end
% 
% return
% 
% function [MDU1,MDU2] = MD(y,M,B,D1,D2,location)
%     DU{1} = reshape(B{1}*D1*flatten(y,1),size(M,1),size(M,2));
%     DU{2} = reshape(B{2}*D2*flatten(y,1),size(M,1),size(M,2));
%     
%     MDU1 = M(:,:,1).*DU{1} + M(:,:,2).*DU{2};
%     MDU2 = M(:,:,3).*DU{1} + M(:,:,4).*DU{2};
% 
%     %MDU1  = interpn(location.ii_v,location.jj_v,MDU1,location.ii_u_extra,location.jj_u_extra,'linear',0);
%     %MDU2  = interpn(location.ii_v,location.jj_v,MDU2,location.ii_u_extra,location.jj_u_extra,'linear',0);
% return