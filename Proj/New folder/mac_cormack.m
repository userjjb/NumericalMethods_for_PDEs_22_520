function Unew = mac_cormack(U,dt);
%
% McCormack's scheme for system of 1D conservation laws
%
% du/dt + d/dx(F + F_A) = 0
% F    = convective flux             
% F_A  = artificial viscous flux    
% 
% The artificial viscosity function have a flag to tell which
% one-sided difference to use.

U = boundary_cond(U);           % Set boundary conditions:

F = flux_function(U);           % Compute flux function:
               
F = F + artificial_visc(-1,U);  % Add artificial viscosity

tmp = dx(1,F);                  % Take a forward difference of the
                                % flux (predictor step)
			       

Up = U - dt*tmp;                % Update solution 
Up = boundary_cond(Up);         % Set boundary conditions on
                                % predictor, Up
F = flux_function(Up);          % Compute flux function using Up
               
F = F + artificial_visc(+1,Up); % Add artificial viscosity

tmp = dx(-1,F);                 % Take a backward difference of the
                                % flux (corrector step)

Unew = 0.5*(U + Up - dt*tmp);   % Update solution 
