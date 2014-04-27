%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Josh Bevan 2013, Updated 2014
%Based on code provided by Prof. Willis, Umass Lowell
%
%Num. Methods for PDEs, 22.520
%Hw 4
 
clc
close all;
clear all;
 
%=========================================================================
% The following parameters are used for defining the mesh of the object
%
% Shape: This selects an ellispoid or a NACA-4 Digit airfoil
%     -- Ellipsoid: Set the shape to a single value (the aspect ratio of
%     the ellipsoid shape. eg: Shape = [1];
%     -- NACA-4-Digit airfoil: Set the shape parameter a 2-valued entry,
%     eg: Shape = [1 0012], where the second value is the 4-digit reference
%     for the NACA airfoil.
%
% DomainSize: This is a basic parameter that defines how far the domain
% should extend (approx.) from your shape. A value of 4 should be
% sufficient, but you can play around with this.
%
% ref: This is the refinement of the mesh. The higher the value, the more
% elements and vertices you will have.
% 
% powerRef: This is the mesh refinment control. The larger the number, the
% more refined the mesh is near the object. A value of 1 produces a uniform
% refiment. A value of 1.75-2.0 should work well for most of your problems.
%
%=========================================================================

                                                    %Shape, DomainSize, ref, powerRef
[TRI, Nodes, Top, Bottom, Left, Right, InnerBound] = getDiscreteGeometry(1, 2, 2, 2);
for i=1:length(InnerBound)
    InnerBoundary(i) = find(and(abs(Nodes(:,1)-InnerBound(i,1))<1e-10,abs(Nodes(:,2)-InnerBound(i,2))<1e-10));
end
NNodes = length(Nodes);
NElem = length(TRI);

%Fancy footwork to put nodal coords in a "nice" format
x=permute(reshape(Nodes(TRI',:)',2,3,numel(TRI)/3),[2,1,3]);
x(:,3,:)=1;
S=diff(x);
TriArea=max(cross(S(1,:,:),S(2,:,:)))/2;

F0 = 0;
F = zeros(NNodes, 1);
for(ie=1:NElem)
    %Calculate gradient coeffs
    Gradient = x(:,:,ie)\eye(3);
    Gradient(3,:) = 0; %Remove unwanted c coeff 
    s(9*ie-8:9*ie) = reshape(TriArea(ie)*(Gradient'*Gradient),9,1);
                      
    F_elemental = (1/3)*TriArea(ie)*[F0 F0 F0]'; %If the f is non-constant this would read in the relevant nodal forces
    F(TRI(ie,:)) = F(TRI(ie,:)) + F_elemental;
end
%Create A matrix, note the vectorized routines for determining global-local
%node mapping for params i and j
A = sparse(...
    reshape(repmat(TRI',3,1),NElem*9,1),...%i
    reshape(repmat(reshape(TRI',1,NElem*3),3,1),NElem*9,1),...%j
    s,NNodes,NNodes,15*NNodes);

% Boundary conditions
Dirichlet = [Left' Right' Top' InnerBoundary];
Neumann = [Bottom];

A(Dirichlet,:)=0;
A((Dirichlet*length(A))+Dirichlet-length(A))=1;
F([Left' Right' Top'])=10;
F([InnerBoundary])=100;

F(Neumann) = 0;
 
% Solve matrix
Sol = A\F;

%_______________POST-PROCESSING__________________
% Gradient calculation
for(ie=1:NElem)
    N1 = TRI(ie,1);
    N2 = TRI(ie,2);
    N3 = TRI(ie,3);
    
    X1 = Nodes(N1,1);
    X2 = Nodes(N2,1);
    X3 = Nodes(N3,1);
    
    Y1 = Nodes(N1,2);
    Y2 = Nodes(N2,2);
    Y3 = Nodes(N3,2);
 
    S1 = Sol(N1);
    S2 = Sol(N2);
    S3 = Sol(N3);
 
    C = [1 X1 Y1; 1 X2 Y2; 1 X3 Y3]\[eye(3)];
    
    Gradient_IE(ie,:) = [(S1*C(2,1)+S2*C(2,2)+S3*C(2,3)), ...
        (S1*C(3,1)+S2*C(3,2)+S3*C(3,3))];
end
 
% Figure 1
figure
trisurf(TRI, Nodes(:,1), Nodes(:,2), Sol-500) % The 500 is random right now
hold on
trisurf(TRI, Nodes(:,1),-Nodes(:,2),Sol-500)
shading interp
title('Scalar Velocity Potential Distribution (Nodal) and the Velocity Vector Field')
view([0 0 1])
axis equal
 
% Figure 2
figure
Vel = ((Gradient_IE(:,1).^2 + Gradient_IE(:,2).^2).^(.5))';
a = trisurf(TRI, Nodes(:,1), Nodes(:,2), Nodes(:,2)*0, Vel);
hold on
set(a,'edgealpha',0)
trisurf(TRI, Nodes(:,1), -Nodes(:,2), -Nodes(:,2)*0, Vel)
title('Velocity Distribution (Centroidal)')
view([0 0 1])
axis equal
