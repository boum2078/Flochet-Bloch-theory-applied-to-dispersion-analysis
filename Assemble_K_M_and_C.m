function [ K, M, C ] = Assemble_K_M_and_C(n_elem_L, L, n_dof, h, t, E, rho, m_discrete , k_discrete, c_discrete )
%ASSAMBLE_K_AND_M Assemble Stiffness, Mass and Damping matrix for unit cells of a 1D
%beam made with a resonator on the last node
%   
%   n_elem_H: number of elements on standing branch
%   H: Heigth of the standing branch
%   n_elem_L: number of elements on the other branches
%   L: Length of the two other branches
%   n_dof: number of dofs for each element-node
%   h: Height of a beam-element
%   t: Thickness of a beam-element
%   theta: angle between L and the vertical, in a counterclockwise direction
%   E: Young modulus [Pa]
%   G: Shear modulus [Pa]
%   rho: Density [kg/m^3]
%   m_discrete: mass of the resonator [kg]
%   k_discrete: stiffness of the resonator [N/m]


% Length of each element
L_elem_L=L/n_elem_L;

% Declaring stiffness, mass and damping matrix; keeping in mind the
% additional dof for the resonators
K=zeros((1+n_elem_L)*n_dof+1, (1+n_elem_L)*n_dof+1);
M=zeros((1+n_elem_L)*n_dof+1, (1+n_elem_L)*n_dof+1);
C=zeros((1+n_elem_L)*n_dof+1, (1+n_elem_L)*n_dof+1);

%  Initialise index counting the number of elements
index=1;

% Computing the mass and stiffness matrix for one element
K_small = Stiffness_EB(L_elem_L,t,h,E);
M_small = Mass_EB(L_elem_L,t,h,rho);

% Pasting the elementary cell n_elem_L times to derive an entire beam
for i=1:n_elem_L,
    K(index:1:(index+2*n_dof-1),index:1:(index+2*n_dof-1))= K(index:1:(index+2*n_dof-1),index:1:(index+2*n_dof-1)) + K_small;
    M(index:1:(index+2*n_dof-1),index:1:(index+2*n_dof-1))= M(index:1:(index+2*n_dof-1),index:1:(index+2*n_dof-1)) + M_small;
    index=index+n_dof;
end

% The resonator mass will be the added to the next node
index_mass = index+n_dof

% Adding the resonator
M(index_mass,index_mass)= M(index_mass,index_mass) + m_discrete;
K([index,index_mass],[index,index_mass])= K([index,index_mass],[index,index_mass])+[k_discrete,-k_discrete;-k_discrete,k_discrete];
C([index,index_mass],[index,index_mass])= C([index,index_mass],[index,index_mass])+[c_discrete,-c_discrete;-c_discrete,c_discrete];
K([index,index_mass],[index,index_mass])


end

