clear;
clf;

n=601;
L=10;
h= L/(n-1);
CFL=0.3
alpha=0.4
t_final=0.004
x=0:h:L;
gamma=1.4
%-------------------------------------------------------------------------%
% Define intial conditions 
p_l= 1e5; %Pressure in left side of shock tube at t=0
p_r= 1e3; %Pressure in right side of shock tube at t=0
rho_l= 1; %Density at left side of shock tube at t=0
rho_r=0.01; %Density at right side of shock tube at t=0
u_l=0; %Velocity in left side of shock tube at t=0
u_r=0; %Velocity in right side of shock tube at t=0
p(1:1:(n+1)/2)=p_l;
p((n+3)/2:1:n)=p_r;
rho(1:1:(n+1)/2)= rho_l;
rho((n+3)/2:1:n)= rho_r;
u(1:1:(n+1)/2)=u_l;
u((N=3)/2:1:n)=u_r;
E= p./((gamma-1)*rho)+0.5*u.^2;%Total Energy
a=sqrt(gamma*p./rho);%Speed of sound
dt=CFL*h/max(a);%Time step based on CFL number
step=0;
%-------------------------------------------------------------------------%
% Time integration begins 
%-------------------------------------------------------------------------%
for t=dt:dt:t_final
    %Define q & F matrix
    q=[rho; rho.*u; rho.*E];
    F=[rho.*u; rho.*u.^2+p; u.*(rho.*E+p)];
    %Calculate q* and flow parameters
    %Calculate F*11/21/07 9:00 AM D:\My Documents\ME 523\proj22\pro_2b.m 2 of 2
    %Calculate artificial viscosity
for i=1:3
 visc(i,1:n-1)=alpha*h^2*rho(1:n-1).*abs((u(2:n)-u(1:n-1))/h).*((u(2:n)-u(1:n-1))
end
%update F* and q matrix
end
%Calculation of flow parameters
%Reference air pressure (N/m^2)
%Reference air density (kg/m^3)
%Entropy w.r.t reference values
%Mass Flow rate per unit area
%-------------------------------------------------------------------------%
% Plot the variables 
%-------------------------------------------------------------------------%
'k' 'X-Coordinate (m)' 'Pressure (Pa)'
'k' 'X-Coordinate (m)' 'Entropy/R gas'
'k' 'X-Coordinate (m)' 'Velocity (m/s)'
'k' 'X-Coordinate (m)' 'Mach number'
'k' 'X-Coordinate (m)' 'Density (kg/m^3)'
'k' 'X-Coordinate (m)' 'Mass Flow (kg/m^2s)'
%---------------------------