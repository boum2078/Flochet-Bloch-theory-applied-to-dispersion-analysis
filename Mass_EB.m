function [ M_global ] = Mass_EB( L,t,h,rho  )

% MASS_EB Computation of the mass matrix of a Euler-Bernouilli
% beam in the global frame
% L         Lengthe of the beam [m]
% h         Height of the beam [m]
% t         Thickness of the beam [m]
% rho       Mass density [kg/m^3]

%% Calculation of derived quantities

Area = t*h;

%% Transformation Matrix from local to world frame
       
T_large=eye(4);         % global coords to local coords transfromation     

%% Euler-Bernouilli Stiffness matrix in Local Frame
             
             
M_local = rho*Area*L*...
                 [156/420,    22*L/420,   	54/420,     -13*L/420;
                  22*L/420,   4*L^2/420,  	13*L/420,   -3*L^2/420;
                  54/420,     13*L/420,   	156/420,    -22*L/420;
                  -13*L/420,  -3*L^2/420, 	-22*L/420,  4*L^2/420];
      
%% Euler-Bernouilli Stiffness matrix in Global Frame

M_global = T_large*M_local*(T_large)';

end
