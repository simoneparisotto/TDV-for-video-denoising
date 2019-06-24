function DNu = discretization(u,params,D1,D2)

%% SIMONE
%%%% first-order %%%%
du = grad(u,D1,D2);

u1 = du(:,:,1);
u2 = du(:,:,2);

%%%% second-order %%%%
du1 = grad(u1,D1,D2);
du2 = grad(u2,D1,D2);

u11 = du1(:,:,1);
u12 = du1(:,:,2);
u21 = du2(:,:,1);
u22 = du2(:,:,2);

%%%% third-order %%%%
du11 = grad(u11,D1,D2);
du12 = grad(u12,D1,D2);
du21 = grad(u21,D1,D2);
du22 = grad(u22,D1,D2);

u111 = du11(:,:,1);
u112 = du11(:,:,2);
u121 = du12(:,:,1);
u122 = du12(:,:,2);
u211 = du21(:,:,1);
u212 = du21(:,:,2);
u221 = du22(:,:,1);
u222 = du22(:,:,2);

switch params.order
    
    case 1
        DNu(:,:,1) = u1;
        DNu(:,:,2) = u2;
    
    case 2
        DNu(:,:,1) = u11;
        DNu(:,:,2) = u12;
        DNu(:,:,3) = u21;
        DNu(:,:,4) = u22;
    case 3
        
        DNu(:,:,1) = u111;
        DNu(:,:,2) = u121;
        DNu(:,:,3) = u211;
        DNu(:,:,4) = u221;
        
        DNu(:,:,5) = u112;
        DNu(:,:,6) = u122;
        DNu(:,:,7) = u212;
        DNu(:,:,8) = u222;
end

end