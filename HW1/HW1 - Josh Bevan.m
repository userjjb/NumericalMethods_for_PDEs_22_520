%Finite difference method for heat transfer in a uniform bar
close all
clear all

%Preallocate
Size=100;
InfNormError = zeros(Size-1,1);

%Define non-changing geometry
LBar = 1;           %Length of bar
Forcing = 1;        %Heat flux on bar

Mode = 2; %1 for normal operation, 2 to examine error

for N=2:Size
if Mode==1
    N=Size;
end
%Geometry
NNodes = N+1;       %number of nodes
DeltaX = LBar/N;    %size of step
X = 0:DeltaX:LBar;  %Nodal locations

%Set up matrices
A = zeros(NNodes, NNodes);
F = zeros(NNodes,1);

%Populate with problem specific equations/values
for(i = 2:(NNodes-1))
    A(i,i) = -2;
    A(i,i-1) = 1;
    A(i,i+1) = 1;
    switch Mode
        case 1
            F(i) = 1;
        case 2
            F(i) = sin(X(i));
    end
end

A = (1./DeltaX^2).*A;

%Boundary conditions
%x=0
A(1, 1) = 1;
F(1) = 0;       %Set T=0
%x=LBar
A(NNodes, NNodes) = 3/(2*DeltaX);       %-
A(NNodes, NNodes-1) = -4/(2*DeltaX);    %Apply back difference for df/dx
A(NNodes, NNodes-2) = 1/(2*DeltaX);     %-
F(NNodes) = 0;  %Set dT/dx=0

%Solve 
Solution = A\F;
switch Mode
    case 1
        Exact = (0.5*(X.^2))-X;
    case 2
        Exact = cos(1)*X - sin(X);
end
Exact = Exact';

if Mode==1
    break
end

%Examine inf-norm error
InfNormError(N-1)=max(abs(Exact-Solution));
end

%Plot results
figure(1)
%plot(X, X*0,'-o')
hold on
plot(X, Solution, '-*r')
temp = plot(X, Exact, 'o');
set(temp,'Color','green')

%Plot error, find slope of loglog curve
if Mode==2
    figure(2)
    N = 2:Size;
    N=N';
    loglog(N,InfNormError, '-r')
    loglogfit = polyfit(log(N),log(InfNormError),1);
    text(10,0.001,num2str(loglogfit(1)))
end

