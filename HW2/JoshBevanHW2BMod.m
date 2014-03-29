clear all
close all
tic

FinalOrder=7; %(2^FinalOrder)+1 number of nodes used along one side
RunMode = 1; %1 for normal operation, 2 to examine error, 3 for execution testing
MapType = 1; %1 for constant transform, 2 for nodal evaluation of transform
SparseType = 2; %1 for indexed population, 2 for apriori creation
if MapType==2 && SparseType==2
    error('Unsupported mode type combination')
end

%Non-loop required initialization
%=========================================================
%Comp domain params
%        xi-eta, xi, eta, offset
CompDomain =   [0 0 0 1;  %a
                0 1 0 1;  %b
                1 1 1 1;  %c
                0 0 1 1]; %d
            
% Physical domain parameters
L = 2;
theta = pi/64;
phi = pi/64;
Forcing = 1;       % Forcing function
%BCs
Left =0;
Right =0;
Top =0;
Bottom =0;

%              a, b           , c                          , d  
PhysDomainX = [0, L*cos(theta), (L*cos(theta))+(L*sin(phi)), L*sin(phi)];
PhysDomainY = [0, L*sin(theta), (L*sin(theta))+(L*cos(phi)), L*cos(phi)];

% Calculate bilinear map coefficients
BiMapX = CompDomain\PhysDomainX';
BiMapY = CompDomain\PhysDomainY';

for Exp=FinalOrder:-1:1
    N=(2^Exp)+1 %Left without semicolon to indicate progress
% Comp domain params
%=========================================================
Nxi = N;  	
Neta = N;  	
NNodes = Nxi*Neta;
F = ones(1,NNodes);
            
Dxi  = 1/(Nxi  - 1);
Deta = 1/(Neta - 1);
            
% Calculate xi, eta and node number
%=========================================================
if SparseType==1
    NodeNumber=zeros(Neta,Nxi);

    counter = 0;
    for j = 1:Neta
        for i = 1:Nxi
            counter = counter + 1;
            NodeNumber(i,j) = counter; 
        end
    end
end

% Implement the mapping 
%=========================================================

%It turns out that the transform coefficients are constant across (i,j) for
%this mapping so we can save a lot of time not calculating them for all 
%i,j) or indexing the matrices in later equations.

% Calculate the transformed equation coefficients
if MapType==1
    a= BiMapX(3)^2+BiMapY(3)^2;
    b= BiMapX(3)*BiMapX(2)+BiMapY(3)*BiMapY(2);
    c= BiMapX(2)^2+BiMapY(2)^2;
    
    alpha = -2*b*BiMapX(1);
    beta  = -2*b*BiMapY(1);

    J= BiMapX(2)*BiMapY(3)-BiMapX(3)*BiMapY(2);

    d= (BiMapY(2)*alpha-BiMapX(2)*beta)/J; 
    e= BiMapX(3)*beta-BiMapY(3)*alpha;
elseif MapType==2
    for j=1:Neta
        for i=1:Nxi
            xi=(i-1)*Neta*Dxi;
            eta=(j-1)*Nxi*Deta;
            a(i,j)= ((BiMapX(1)*xi+BiMapX(3))^2)+((BiMapY(1)*xi+BiMapY(3))^2);
            b(i,j)= ((BiMapX(1)*xi+BiMapX(3))*(BiMapX(1)*eta+BiMapX(2)))+((BiMapY(1)*xi+BiMapY(3))*(BiMapY(1)*eta+BiMapY(2)));
            c(i,j)= ((BiMapX(1)*eta+BiMapX(2))^2)+((BiMapY(1)*eta+BiMapY(2))^2);

            alpha = -2*b(i,j)*BiMapX(1);
            beta  = -2*b(i,j)*BiMapY(1);

            J(i,j)= ((BiMapX(1)*eta+BiMapX(2))*(BiMapY(1)*xi+BiMapY(3)))-((BiMapX(1)*xi+BiMapX(3))*(BiMapY(1)*eta+BiMapY(2)));

            d(i,j)= (((BiMapY(1)*eta+BiMapY(2))*alpha)-((BiMapX(1)*eta+BiMapX(2))*beta))/J(i,j); 
            e(i,j)= ((BiMapX(1)*xi+BiMapX(3))*beta)-((BiMapY(1)*xi+BiMapY(3))*alpha);
       end
    end
end

% Set-up finite difference matrix 
%=========================================================
if MapType==1
    corner =(2*b)/(4*Deta*Dxi)*(-1/J^2);
    ip =(a/Dxi^2 +e/(2*Dxi))*(-1/J^2);
    in =(a/Dxi^2 -e/(2*Dxi))*(-1/J^2);
    jp =(c/Deta^2 +d/(2*Deta))*(-1/J^2);
    jn =(c/Deta^2 -d/(2*Deta))*(-1/J^2);
    center =((-2*c)/Deta^2 -(2*a)/Dxi^2)*(-1/J^2);

    if SparseType==1
        A = spalloc(NNodes,NNodes, NNodes*10);	% allocate the matrix
        for j = 2:(Neta - 1)
            for i = 2:(Nxi-1)
                % Populate A-Matrix
                %================================================
                NodeN = NodeNumber(i,j);
                A(NodeN, NodeNumber(i+0,j+0)) =  center; 

                A(NodeN, NodeNumber(i+1,j+0)) =  ip;
                A(NodeN, NodeNumber(i-1,j+0)) =  in;

                A(NodeN, NodeNumber(i+0,j+1)) =  jp;
                A(NodeN, NodeNumber(i+0,j-1)) =  jn;

                A(NodeN, NodeNumber(i+1,j+1)) =  -corner;
                A(NodeN, NodeNumber(i+1,j-1)) =  corner;
                A(NodeN, NodeNumber(i-1,j+1)) =  corner;
                A(NodeN, NodeNumber(i-1,j-1)) =  -corner;
            end
        end
    elseif SparseType==2
        %Create matrix containing diagonals to be put in A matrix
        %Brother are directly adjacent diags(nodes along i in stencil), cousins are the
        %grouped diags(each level along j in stencil)
        %"Magic" numbers are static and based on a 3x3 stencil
        DiagValue = ones(NNodes,9);
        for cousin =0:2
            for brother =0:2
                for iter =1:Neta-3
                    R =Nxi*iter-1:Nxi*iter;
                    %Snips out gaps amongst diags in A
                    DiagValue(R+brother+cousin*Nxi,cousin*3+brother+1) =0;
                end
                %Cuts off top and bottom of A
                DiagValue([1:brother+cousin*Nxi, NNodes-(2-cousin)*Nxi+brother-1:NNodes],cousin*3+brother+1) =0;
            end
        end
        %Set diags to correct static value
        DiagValue(:,1) =DiagValue(:,1).*-corner;
        DiagValue(:,2) =DiagValue(:,2).*jn;
        DiagValue(:,3) =DiagValue(:,3).*corner;
        DiagValue(:,4) =DiagValue(:,4).*in;
        DiagValue(:,5) =DiagValue(:,5).*center;
        DiagValue(:,6) =DiagValue(:,6).*ip;
        DiagValue(:,7) =DiagValue(:,7).*corner;
        DiagValue(:,8) =DiagValue(:,8).*jp;
        DiagValue(:,9) =DiagValue(:,9).*-corner;
        %Set-up BCs
            %BC finite differences
        DiagValue(1:Nxi,5) = 1; %Bottom
        DiagValue(NNodes-Nxi+1:NNodes,5) = 1; %Top
        DiagValue(Nxi+1:Nxi:NNodes-Nxi+1,5) = 1; %Left
        DiagValue(2*Nxi:Nxi:NNodes-Nxi,5) = 1; %Right
            %BC values
        F(1:Nxi)=Bottom;
        F(NNodes-Nxi+1:NNodes) =Top;
        F(Nxi+1:Nxi:NNodes-Nxi+1) =Left;
        F(2*Nxi:Nxi:NNodes-Nxi) =Right;
                
        %Create sparse matrix with diags defined in DiagValue along diag numbers
        %defined in DiagLoc
        %Based on 3x3 stencil
        DiagLoc =[-Nxi-1 -Nxi -Nxi+1 -1 0 1 Nxi-1 Nxi Nxi+1];
        A=spdiags(DiagValue, DiagLoc, NNodes, NNodes);
    end
elseif MapType==2
    A = spalloc(NNodes,NNodes, NNodes*10);	% allocate the matrix
    for j = 2:(Neta - 1)
        for i = 2:(Nxi-1)
            % Populate A-Matrix
            %================================================
            NodeN = NodeNumber(i,j);
            A(NodeN, NodeNumber(i+0,j+0)) =  ((-2*c(i,j))/Deta^2 -(2*a(i,j))/Dxi^2)     *(-1./J(i,j)^2);   

            A(NodeN, NodeNumber(i+1,j+0)) =  (a(i,j)/Dxi^2 +e(i,j)/(2*Dxi))             *(-1./J(i,j)^2);
            A(NodeN, NodeNumber(i-1,j+0)) =  (a(i,j)/Dxi^2 -e(i,j)/(2*Dxi))             *(-1./J(i,j)^2);

            A(NodeN, NodeNumber(i+0,j+1)) =  (c(i,j)/Deta^2 +d(i,j)/(2*Deta))           *(-1./J(i,j)^2);
            A(NodeN, NodeNumber(i+0,j-1)) =  (c(i,j)/Deta^2 -d(i,j)/(2*Deta))           *(-1./J(i,j)^2);

            A(NodeN, NodeNumber(i+1,j+1)) =  ((-2*b(i,j))/(4*Deta*Dxi))                 *(-1./J(i,j)^2);
            A(NodeN, NodeNumber(i+1,j-1)) =  ((2*b(i,j))/(4*Deta*Dxi))                  *(-1./J(i,j)^2);
            A(NodeN, NodeNumber(i-1,j+1)) =  ((2*b(i,j))/(4*Deta*Dxi))                  *(-1./J(i,j)^2);
            A(NodeN, NodeNumber(i-1,j-1)) =  ((-2*b(i,j))/(4*Deta*Dxi))                 *(-1./J(i,j)^2);
        end
    end
end

%Populate forcing vector
%================================================
F=F.*Forcing;

% Implement Boundary conditions
%=========================================================
if SparseType==1
    for(i = 1:(Nxi))
        j = 1;
        ij = NodeNumber(i,j);
        A(ij,ij) = Bottom;
        F(ij) = 0;

        j = Neta;
        ij = NodeNumber(i,j);
        A(ij,ij) = Top;
        F(ij) = 0;
    end

    for(j = 2:(Neta-1))
        i = 1;
        ij = NodeNumber(i,j);
        A(ij,ij) = Left;
        F(ij) = 0;

        i = Nxi;
        ij = NodeNumber(i,j);
        A(ij,ij) = Right;
        F(ij) = 0;
    end
end

% Solve the linear system
%=========================================================
Solution = A\F';
MaxDeflect = max(max(Solution));
if Exp==FinalOrder
    Exact=MaxDeflect;
else
    InfNormError(Exp)=abs(Exact-MaxDeflect);
end
Saved(Exp)=MaxDeflect; %Look for spurious behavior, this should stay mostly constant

if RunMode~=2 %Break after 1 loop for normal mode
    break
end
end

%Plot solved function for a particular N
if RunMode==1
    Deflection = reshape(Solution, Nxi, Neta);
    figure
    surf(0:Dxi*Neta:Nxi,0:Deta*Neta:Neta,Deflection)
    colorbar
    axis equal
    view([0 0 1])
    title('Computational Domain')
    
    %Find mapped locations to plot physical domain
    for j=1:Neta
        for i=1:Nxi
            xi=(i-1)*Dxi;
            eta=(j-1)*Deta;
            X(i,j) = (BiMapX(1)*(xi*eta))+ (BiMapX(2)*xi) + (BiMapX(3)*eta) + BiMapX(4);
            Y(i,j) = (BiMapY(1)*(xi*eta))+ (BiMapY(2)*xi) + (BiMapY(3)*eta) + BiMapY(4);
        end
    end
    
    figure
    surf(X,Y,Deflection)
    colorbar
    axis equal
    view([0 0 1])
    title('Physical Domain')
end

%Plot error, find slope of loglog curve
if RunMode==2
    N = 1:FinalOrder-1; %Exclude final order itself
    N = (2.^N)+1;
    loglog(N',InfNormError, '-r')
    loglogfit = polyfit(log(N),log(InfNormError),1);
    text(N(floor(FinalOrder/2)),InfNormError(floor(FinalOrder/2)),num2str(loglogfit(1)))
end

toc
beep



