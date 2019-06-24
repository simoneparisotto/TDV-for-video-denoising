clear
close all
clc

addpath ./../lib/lib_TDV

addpath ./../lib/export_fig-master/

% plot parameters
lw   = 3;
aw   = 0.3;
fsz  = 25;
fsza = 18;

%% ROTATION MATRICES (the number defines the rotation axis 1=i, 2=j, 3=k)
R1 = @(theta) [1 0 0                    ;  0 cos(theta) -sin(theta) ; 0 sin(theta) cos(theta)];
R2 = @(theta) [cos(theta) 0 -sin(theta) ;  0 1 0                    ; sin(theta) 0 cos(theta)];
R3 = @(theta) [cos(theta) -sin(theta) 0 ; sin(theta) cos(theta) 0   ; 0 0 1                  ];

%% PROJECTION MATRICES
P12 = [1 0 0; 0 1 0; 0 0 0];
P13 = [1 0 0; 0 0 0; 0 0 1];
P23 = [0 0 0; 0 1 0; 0 0 1];

%% CARTESIAN AXES [x,y,t]
xyt = eye(3);

%% ORTHOGONAL VECTORS [u,v,z]: eigenvectors + eigenvalues
% define the main unitary direction wrt the cartesian axes
u = [-0.5 0.8 2];
u=u(:).'/norm(u);
vz=null(u).';

% eigenvectors
uvz = [u;vz];

% possible eigenvalues
lambda = [1 0.5 0.1];

%% CREATE 3D TENSOR FOR ELLIPSOID
S123 = lambda(1)*kron(uvz(:,1),uvz(:,1).') + lambda(2)*kron(uvz(:,2),uvz(:,2).') + lambda(3)*kron(uvz(:,3),uvz(:,3).');
[uvz,lambda] = eig(S123,'vector');

% DEFINE THE ELLIPSOID described by the eigenvectors uvx and eigenvalues lambda
uvz = uvz * diag(sqrt(lambda)); 

 % reconstruction for test purposes
S123bis = kron(uvz(:,1),uvz(:,1).') + kron(uvz(:,2),uvz(:,2).') + kron(uvz(:,3),uvz(:,3).');

%% PROJECT ELLIPSE AXES uvz ON xyt
P12x = P12*uvz;
P13x = P13*uvz;
P23x = P23*uvz;

%% FIND ROTATION MATRICES
% rotate {x,y} around axes t
a      = uvz(:,3); 
b      = [0 0 1];
na     = a./norm(a);
nb     = b./norm(b);
alpha3 = rad2deg(atan2(norm(cross(a, b)), dot(a, b)));
vv     = cross(na, nb);
skew   = [0, -vv(3), vv(2); vv(3), 0, -vv(1); -vv(2), vv(1), 0];
RR3    = eye(3) + skew + skew ^ 2 * (1 - dot(na, nb)) / (norm(vv))^2;
R3x    = RR3 * uvz;

% rotate {x,t} around axes y
a      = R3x(:,1);
b      = [0 1 0];
na     = a./norm(a);
nb     = b./norm(b);
alpha2 = rad2deg(atan2(norm(cross(a, b)), dot(a, b)));
vv     = cross(na, nb);
skew   = [0, -vv(3), vv(2); vv(3), 0, -vv(1); -vv(2), vv(1), 0];
RR2    = eye(3) + skew + skew ^ 2 * (1 - dot(na, nb)) / (norm(vv))^2;
R2R3x  = RR2 * R3x;

RR    = RR2*RR3;
invRR = RR3'*RR2';

% ROTATE uvz' BACK IN THE ORIGINAL PLACE
iR3R2x = invRR*R2R3x;

%% COMPUTE ELLIPSOID
N = 80; % number of points in the grid
[ELL(:,:,1),ELL(:,:,2),ELL(:,:,3)] = ellipsoid(0,0,0,sqrt(lambda(2)),sqrt(lambda(1)),sqrt(lambda(3)),N);
ell = reshape(ELL,[],3);

% rotate ellipse
rotell = (invRR*ell.').';

% project ellipse onto cartesian axes
P12rotell = (P12*rotell.').';
P13rotell = (P13*rotell.').';
P23rotell = (P23*rotell.').';

ROTELL(:,:,1)    = reshape(rotell(:,1),N+1,N+1);
ROTELL(:,:,2)    = reshape(rotell(:,2),N+1,N+1);
ROTELL(:,:,3)    = reshape(rotell(:,3),N+1,N+1);
P12ROTELL(:,:,1) = reshape(P12rotell(:,1),N+1,N+1);
P12ROTELL(:,:,2) = reshape(P12rotell(:,2),N+1,N+1);
P12ROTELL(:,:,3) = reshape(P12rotell(:,3),N+1,N+1);
P13ROTELL(:,:,1) = reshape(P13rotell(:,1),N+1,N+1);
P13ROTELL(:,:,2) = reshape(P13rotell(:,2),N+1,N+1);
P13ROTELL(:,:,3) = reshape(P13rotell(:,3),N+1,N+1);
P23ROTELL(:,:,1) = reshape(P23rotell(:,1),N+1,N+1);
P23ROTELL(:,:,2) = reshape(P23rotell(:,2),N+1,N+1);
P23ROTELL(:,:,3) = reshape(P23rotell(:,3),N+1,N+1);

%% DEFINE THE TENSORS OF THE PROJECTED ELLIPSES AND EIGENDECOPOSITION

% COMPUTE THE TENSOR OF THE ELLIPSOID PROJECTED ONTO {x,y}
S12 = zeros(3);
% work with boundary points only
K = convhulln(P12rotell(:,[1,2]));  
K = unique(K(:));  
Q = P12rotell(K,1:2).';
[S12([1,2],[1,2]),c] = MinVolEllipse(Q, 0.001);
[e12,l12] = eig(S12,'vector'); 
e12 = e12 * diag(sqrt(1./l12));

e12n = e12; 
e12n(:,1) = e12n(:,1)./sqrt(1./l12(2));
e12n(:,2) = e12n(:,2)./sqrt(1./l12(2));

%% FIGURES
if true
    figure('Position',[1 1 1000 1000])
    subplot(3,3,1)
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],uvz(1,:),uvz(2,:),uvz(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectors
    quiver3(0,0,0,u(1),u(2),u(3),'linewidth',lw,'MaxHeadSize',aw,'color','m','AutoScale','off') % main direction (associated to largest eigenvalues)
    title('Eigenvectors (eigbasis)')
    axis square
    grid on
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    view([130 20])
    
    % project on (x,y)
    subplot(3,3,4)
    patch( [1 -1 -1 1] , [-1 -1 1 1], [0 0 0 0], 'red','FaceAlpha',.2) % (x,y)
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    quiver3([0,0,0],[0,0,0],[0,0,0],uvz(1,:),uvz(2,:),uvz(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectors
    quiver3(0,0,0,u(1),u(2),u(3),'linewidth',lw,'MaxHeadSize',aw,'color','m','AutoScale','off') % main direction (associated to largest eigenvalues)
    quiver3([0,0,0],[0,0,0],[0,0,0],P12x(1,1:3),P12x(2,1:3),P12x(3,1:3),'linewidth',lw,'MaxHeadSize',aw,'color','g','AutoScale','off')
    axis square
    grid on
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    title('Project eigbasis on (x,y)')
    view([130 20])
    
    % project on (x,z)
    subplot(3,3,5)
    patch( [1 -1 -1 1] , [0 0 0 0], [1 1 -1 -1], 'red','FaceAlpha',.2) % (x,z)
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    quiver3([0,0,0],[0,0,0],[0,0,0],uvz(1,:),uvz(2,:),uvz(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectors
    quiver3(0,0,0,u(1),u(2),u(3),'linewidth',lw,'MaxHeadSize',aw,'color','m','AutoScale','off') % main direction (associated to largest eigenvalues)
    quiver3([0,0,0],[0,0,0],[0,0,0],P13x(1,:),P13x(2,:),P13x(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','g','AutoScale','off')
    axis square
    grid on
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    title('Project eigbasis on (x,z)')
    view([130 20])
    
    % project on (y,z)
    subplot(3,3,6)
    patch( [0 0 0 0], [1 -1 -1 1] , [-1 -1 1 1], 'red','FaceAlpha',.2) % (y,z) % (y,z)
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    quiver3([0,0,0],[0,0,0],[0,0,0],uvz(1,:),uvz(2,:),uvz(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectors
    quiver3(0,0,0,u(1),u(2),u(3),'linewidth',lw,'MaxHeadSize',aw,'color','m','AutoScale','off') % main direction (associated to largest eigenvalues)
    quiver3([0,0,0],[0,0,0],[0,0,0],P23x(1,:),P23x(2,:),P23x(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','g','AutoScale','off')
    axis square
    grid on
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    title('Project eigbasis on (y,z)')
    view([130 20])
    
    subplot(3,3,7)
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],uvz(1,:),uvz(2,:),uvz(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectors
    axis square
    grid on
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    title('Initial status')
    view([130 20])
    
    subplot(3,3,8)
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],R2R3x(1,:),R2R3x(2,:),R2R3x(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectorstitle('Eigenvectors (eigbasis)')
    axis square
    grid on
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    title('Rotate the uvzrot = RR2*RR3*uvz')
    view([130 20])
    
    subplot(3,3,9)
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],iR3R2x(1,:),iR3R2x(2,:),iR3R2x(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectorstitle('Eigenvectors (eigbasis)')
    axis square
    grid on
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    title('Rotate the inv(RR3)*inv(RR2)*uvzrot')
    view([130 20])
end

%% PLOT ELLIPSOIDS
if true
    figure('Position',[1 1 1000 1000])
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],uvz(1,:),uvz(2,:),uvz(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectors
    quiver3([0,0,0],[0,0,0],[0,0,0],R2R3x(1,:),R2R3x(2,:),R2R3x(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','c','AutoScale','off') % orthonormal eigenvectors
    surf(ELL(:,:,1),ELL(:,:,2),ELL(:,:,3),'Facecolor','cyan','FaceAlpha',0.3,'FaceLighting','gouraud')
    shading interp
    camlight
    axis square
    grid on
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    title('Initial status')
    view([130 20])
    
    h2 = figure('Position',[1 1 1000 1000]);
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],iR3R2x(1,:),iR3R2x(2,:),iR3R2x(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectors
    surf(ROTELL(:,:,1),ROTELL(:,:,2),ROTELL(:,:,3),'Facecolor','cyan','FaceAlpha',0.3,'FaceLighting','gouraud')
    colormap([0 0 1])
    shading interp
    camlight
    axis square
    grid on
    set(gca,'fontsize',fsza)
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    view([130 20])
    set(gcf, 'Color', 'w');
    export_fig(h2,'ellipsoid_rotated.png')
    title('Rotate back ellipsoid')
    
    h3 = figure('Position',[1 1 1000 1000]);
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],iR3R2x(1,:),iR3R2x(2,:),iR3R2x(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectors
    quiver3([0,0,0],[0,0,0],[0,0,0],e12(1,:),e12(2,:),e12(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','g','AutoScale','off') % orthonormal eigenvectors
    surf(P12ROTELL(:,:,1),P12ROTELL(:,:,2),P12ROTELL(:,:,3),'Facecolor','red','FaceAlpha',0.3)
    shading interp
    surf(P13ROTELL(:,:,1),P13ROTELL(:,:,2),P13ROTELL(:,:,3),'Facecolor','green','FaceAlpha',0.3)
    shading interp
    surf(P23ROTELL(:,:,1),P23ROTELL(:,:,2),P23ROTELL(:,:,3),'Facecolor','blue','FaceAlpha',0.3)
    colormap([1 0 0])
    shading interp
    axis square
    grid on
    set(gca,'fontsize',fsza)
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    view([130 20])
    set(gcf, 'Color', 'w');
    export_fig(h3,'ellipsoid_rotated_projected.png')
    title('Projected ellipsoid')
    
    
    h4 = figure('Position',[1 1 1000 1000]);
    quiver3([0,0,0],[0,0,0],[0,0,0],xyt(1,:),xyt(2,:),xyt(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','r','AutoScale','off') % cartesian axes
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],iR3R2x(1,:),iR3R2x(2,:),iR3R2x(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','b','AutoScale','off') % orthonormal eigenvectors
    quiver3([0,0,0],[0,0,0],[0,0,0],e12n(1,:),e12n(2,:) ,e12n(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','m','AutoScale','off') % orthonormal eigenvectors
    quiver3([0,0,0],[0,0,0],[0,0,0],e12(1,:),e12(2,:),e12(3,:),'linewidth',lw,'MaxHeadSize',aw,'color','g','AutoScale','off') % orthonormal eigenvectors
    surf(P12ROTELL(:,:,1),P12ROTELL(:,:,2),P12ROTELL(:,:,3),'Facecolor','red','FaceAlpha',0.3)
    shading interp
    surf(P13ROTELL(:,:,1),P13ROTELL(:,:,2),P13ROTELL(:,:,3),'Facecolor','green','FaceAlpha',0.3)
    shading interp
    surf(P23ROTELL(:,:,1),P23ROTELL(:,:,2),P23ROTELL(:,:,3),'Facecolor','blue','FaceAlpha',0.3)
    colormap([1 0 0])
    shading interp
    axis square
    grid on
    set(gca,'fontsize',fsza)
    xlim([-1,1]), xlabel('x','Fontsize',fsz)
    ylim([-1,1]), ylabel('y','Fontsize',fsz)
    zlim([-1,1]), zlabel('t','Fontsize',fsz)
    view([0 90])
    set(gcf, 'Color', 'w');
    export_fig(h4,'ellipsoid_rotated_projected_top.png')
    title('Projected ellipsoid')
end