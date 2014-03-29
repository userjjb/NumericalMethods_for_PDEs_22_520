%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Josh Bevan 2013
%Based on code by Prof. WIllis, Umass Lowell
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
% refiment.
%
%=========================================================================

Shape = [1 0045];
DomainSize = 6;
powerRef = 2;
ref=3;

[TRI, Nodes, Top, Bottom, Left, Right, InnerBoundary] = getDiscreteGeometry(Shape, DomainSize, ref, powerRef);
%===
% plot the triangulation
%===
trimesh(TRI, Nodes(:,1),Nodes(:,2))
axis equal
pause(0.1)

%===
% Determine the number of nodes and elements in the domain
%===
NNodes = length(Nodes)
NElem = length(TRI);

%===
% Initialize the A-matrix and RHS to zero valued containers (could use
% spalloc)
%===
A = spalloc(NNodes, NNodes, 12*NNodes);
F = zeros(NNodes, 1);
F0 = 0;

%===
% cycle through the elements to build the A-matrix and f-vector
%===
for(ie=1:NElem)
    element_number = ie;

    N1 = TRI(ie,1);
    N2 = TRI(ie,2);
    N3 = TRI(ie,3);
    
    X1 = Nodes(N1,1);
    X2 = Nodes(N2,1);
    X3 = Nodes(N3,1);
    
    Y1 = Nodes(N1,2);
    Y2 = Nodes(N2,2);
    Y3 = Nodes(N3,2);
    
    TriArea = calculateArea(X1, X2, X3, Y1, Y2, Y3);
    
    %Calculate gradient coeffs
    Plane = [X1 Y1 1;...
             X2 Y2 1;...
             X3 Y3 1];
    Gradient(:,1) = Plane\[1;0;0];
    Gradient(:,2) = Plane\[0;1;0]; 
    Gradient(:,3) = Plane\[0;0;1];
    Gradient(3,:) = 0; %Remove unwanted c coeff
                        
    gN1gN1 = dot(Gradient(:,1),Gradient(:,1));
    gN2gN2 = dot(Gradient(:,2),Gradient(:,2));
    gN3gN3 = dot(Gradient(:,3),Gradient(:,3));

    gN1gN2 = dot(Gradient(:,1),Gradient(:,2));
    gN1gN3 = dot(Gradient(:,1),Gradient(:,3));            
    gN2gN3 = dot(Gradient(:,2),Gradient(:,3));

    A_elemental = TriArea*[gN1gN1  gN1gN2  gN1gN3;...
                           gN1gN2  gN2gN2  gN2gN3;...
                           gN1gN3  gN2gN3  gN3gN3];
               
              
    F_elemental = (1/3)*TriArea*[F0 F0 F0]; %If the f is non-constant this would read in the relevant nodal forces
 
    %Note that manually poking values to A and F are slow, it would be better to
    %pre-generate the diagonals and create the sparse matrix with spdiags()  
    A(N1, N1) = A(N1, N1) + A_elemental(1,1);
    A(N2, N1) = A(N2, N1) + A_elemental(2,1);
    A(N3, N1) = A(N3, N1) + A_elemental(3,1);

    A(N1, N2) = A(N1, N2) + A_elemental(1,2);
    A(N2, N2) = A(N2, N2) + A_elemental(2,2);
    A(N3, N2) = A(N3, N2) + A_elemental(3,2);
    
    A(N1, N3) = A(N1, N3) + A_elemental(1,3);
    A(N2, N3) = A(N2, N3) + A_elemental(2,3);
    A(N3, N3) = A(N3, N3) + A_elemental(3,3);

    F(N1) = F(N1) + F_elemental(1);
    F(N2) = F(N2) + F_elemental(2);
    F(N3) = F(N3) + F_elemental(3);
end

%===
% Boundary conditions
%===
%Neumann BCs on Top/Bottom/Airfoil
for i=1:length(Top)
    F(Top(i)) = 0;
end
for i=1:length(Bottom)
    F(Bottom(i)) = 0;
end
%Airfoil boundary
for i=1:length(InnerBoundary)
    searcher = InnerBoundary(i,1);
    found = find(abs(Nodes(:,1)-searcher)<1e-10);
    F(found) = 0;
end

%Dirichlet BCs on Left/Right
for i=1:length(Left)
    %Left
    A(Left(i),:) = 0;
    A(Left(i),Left(i)) = 1;
    F(Left(i)) = Nodes(Left(i),1);
    %Right
    A(Right(i),:) = 0;
    A(Right(i),Right(i)) = 1;
    F(Right(i)) = Nodes(Right(i),1);
end

%===
% Solve the problem
%===
Sol = A\F;

%===
% Gradient calculation
%===
for(ie=1:NElem)
    element_number = ie;

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
    
    Centroid(ie,:) = [(X1+X2+X3)/3, (Y1+Y2+Y3)/3];
end

% ===
% Figure 1
% ===
figure
trisurf(TRI, Nodes(:,1), Nodes(:,2), Sol)
hold on
trisurf(TRI, Nodes(:,1),-Nodes(:,2),Sol)
quiver(Centroid(:,1), Centroid(:,2), Gradient_IE(:,1), Gradient_IE(:,2),.75,'k')
quiver(Centroid(:,1), -Centroid(:,2), Gradient_IE(:,1), -Gradient_IE(:,2),.75,'k')
shading interp
title('Scalar Velocity Potential Distribution (Nodal) and the Velocity Vector Field')
view([0 0 1])
axis equal


% ===
% Figure 2
% ===
figure
Vel = ((Gradient_IE(:,1).^2 + Gradient_IE(:,2).^2).^(.5))';
a = trisurf(TRI, Nodes(:,1), Nodes(:,2), Nodes(:,2)*0, Vel);
hold on
set(a,'edgealpha',0)
trisurf(TRI, Nodes(:,1), -Nodes(:,2), -Nodes(:,2)*0, Vel)
title('Velocity Distribution (Centroidal)')
view([0 0 1])
axis equal