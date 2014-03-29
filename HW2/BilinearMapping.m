clear all
close all

L = 2;
theta = pi/3.5;
phi = pi/12;
%        xi-eta, xi, eta, offset
CompDomain =   [0 0 0 1;  %a
                0 1 0 1;  %b
                1 1 1 1;  %c
                0 0 1 1]; %d
            
%              a, b           , c                          , d  
PhysDomainX = [1, L*cos(theta)+1, (L*cos(theta))+(L*sin(phi))+1, L*sin(phi)+1];
PhysDomainY = [0, L*sin(theta), (L*sin(theta))+(L*cos(phi)), L*cos(phi)];
            
BilinearMapX = CompDomain\PhysDomainX';
BilinearMapY = CompDomain\PhysDomainY';

Ntests = 10;

axis equal
hold on
TestPointsXi = rand(Ntests, 1)*max(CompDomain(:,2))+min(CompDomain(:,2));
TestPointsEta = rand(Ntests, 1)*max(CompDomain(:,3))+min(CompDomain(:,2));
XiEta = TestPointsXi.*TestPointsEta;
for iter=1:Ntests
    MappedX(iter) = (BilinearMapX(1)*XiEta(iter))+ (BilinearMapX(2)*TestPointsXi(iter))+ (BilinearMapX(3)*TestPointsEta(iter))+ +BilinearMapX(4);
    MappedY(iter) = (BilinearMapY(1)*XiEta(iter))+ (BilinearMapY(2)*TestPointsXi(iter))+ (BilinearMapY(3)*TestPointsEta(iter))+ +BilinearMapY(4);
end

axis equal
fill(PhysDomainX,PhysDomainY,[0 .6 1])
fill(CompDomain(:,2),CompDomain(:,3),[.7 .2 .2])

plot(TestPointsXi,TestPointsEta,'ro')
plot(MappedX,MappedY,'bo')