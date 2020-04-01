function [ K_global ] = Stiffness_EB(L,t,h,E)

% STIFFNESS_EB Computation of the stiffness matrix of a Euler-Bernouilli
% beam in the global frame
% L         Lengthe of the beam [m]
% h         Height of the beam [m]
% t         Thickness of the beam [m]
% E         Young modulus of the beam [Pa]
% G         Shear modulus of the beam [Pa]
% angle     Angle between the beam and the globel x-axis [rad] (counterclockwise) 

%% Calculation of derived quantities

% I_zz = pi*r^4/4;

I_zz = t*h^3/12;
if h>t
    J=h*t^3*(1/3-0.21*t/h*(1-t^4/(12*h^4)));
elseif t>=h
    J=t*h^3*(1/3-0.21*h/t*(1-h^4/(12*t^4)));
else
    
end

%% Transformation Matrix from local to world frame

T_large=eye(4);       % global coords to local coords transfromation

%% Euler-Bernouilli Stiffness matrix in Local Frame

             
K_local = 1/L^3*[12*E*I_zz,  6*E*I_zz*L,     -12*E*I_zz,     6*E*I_zz*L;
                 
                 6*E*I_zz*L, 4*E*I_zz*L^2,   -6*E*I_zz*L,    2*E*I_zz*L^2;
                 -12*E*I_zz, -6*E*I_zz*L,    12*E*I_zz,      -6*E*I_zz*L;
                
                 6*E*I_zz*L, 2*E*I_zz*L^2,   -6*E*I_zz*L,    4*E*I_zz*L^2];
             
             
                          
%% Euler-Bernouilli Stiffness matrix in Global Frame

K_global = T_large*K_local*(T_large)';

end

