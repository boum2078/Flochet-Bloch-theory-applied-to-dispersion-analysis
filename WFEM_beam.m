clear all
close all
clc
%% Definition of the bar

% Length of the unit cell
L=0.1;                   
% height of the beam
h=L/10;
%Thickness of the beam
t=h;

% Young's modulus [Pa]
E = 2.1e11;
% Density [kg/m^3]
rho = 7800;             

%% Definition of the added disturbance

% Added mass to the beam
m_ratio=0.2;
% Resonator = true if resonator, if false = point mass
resonator=1;
% Frequency of the resonator
f_band_low_dimensionless=0.5001;
% Percent proportional damping in the resonator
ksi=0;

%% Derived variables
% Second moment of Inertia
I_zz = t*h^3/12;
% Surface area of the beam
Area = h*t;
% Wave speed in the beam
c=sqrt(E/rho);
% Massa of the original beam
m_beam=rho*Area*L;
% Natural frequency of the beam (standing wave across 1 cell)
nat_frequency=sqrt(c^2.*(pi/L).^4*I_zz/Area)/(2*pi);
% Mass of the point mass
m_discrete=m_ratio*m_beam;
% Stiffnes of the resonator
if resonator
    k_discrete=(f_band_low_dimensionless*nat_frequency*2*pi)^2*m_discrete;
else
    k_discrete=0;
end
% Damping value
c_discrete=ksi*2*sqrt(m_discrete*k_discrete);

%% Building the model

% Frequency ratio upto which the calculations are performed
max_f_ratio=1.5;

% Number of elements for the beam
n_elem_L = 10;

% Number of dofs per node
n_dof=2;

% Build the model
% this part can be replaced by any FE model of a unit cell you would like to extract its wave based properties.
% such models can be meshed in a third party software and expported to Matlab
[ K, M, C ] = Assemble_K_M_and_C(n_elem_L, L, n_dof, h, t, E, rho, m_discrete , k_discrete, c_discrete );

%% Solving system

% Setting the radial frequencies for which the EVP will be solved
omega=linspace(0,2000,500)*2*pi;

% left side of UC, the dofs go from 1 to n_ql
n_ql=n_dof;
% right side of UC, the dofs go from n_qr to 1
n_qr=n_dof;
if resonator
    figure,spy(M)
    % shift the resonator location so that the last elements are boundary
    % elements
    K=K([1:end-(n_qr+1),end,(end-n_qr):(end-1)],[1:end-(n_qr+1),end,(end-n_qr):(end-1)]);
    M=M([1:end-(n_qr+1),end,(end-n_qr):(end-1)],[1:end-(n_qr+1),end,(end-n_qr):(end-1)]);
    C=C([1:end-(n_qr+1),end,(end-n_qr):(end-1)],[1:end-(n_qr+1),end,(end-n_qr):(end-1)]);
%     figure,spy(M)
else
    % if no resonator, delete the resonator dof and put the mass on the
    % last element
    K=K(1:end-1,1:end-1);
    M=M(1:end-1,1:end-1);
    C=C(1:end-1,1:end-1);
    M(end-1,end-1)=M(end-1,end-1)+m_discrete;
end

% Declaring the solution vector 
k_times_L=zeros(n_dof*2,length(omega));

Nodes.Left=[1:n_ql];
Nodes.Right=[length(K)-n_qr+1:length(K)];
Nodes.Intern=setdiff(1:length(K),[Nodes.Left,Nodes.Right]);

for ind=1:length(omega)
    % Calculating the dynamic stiffness matrix
    D=K+sqrt(-1)*omega(ind)*C-omega(ind)^2*M;
    
    % Orderint the matrix according to the node set a submatrix belongs to
    D_LL_t=D(Nodes.Left,Nodes.Left);
    D_IL_t=D(Nodes.Intern,Nodes.Left);
    D_RL_t=D(Nodes.Right,Nodes.Left);
    
    D_LI_t=D(Nodes.Left,Nodes.Intern);
    D_II_t=D(Nodes.Intern,Nodes.Intern);
    D_RI_t=D(Nodes.Right,Nodes.Intern);
    
    D_LR_t=D(Nodes.Left,Nodes.Right);
    D_IR_t=D(Nodes.Intern,Nodes.Right);
    D_RR_t=D(Nodes.Right,Nodes.Right);
    
    
    
    % Condensation of the internal DOFs 
    D_LL=D_LL_t-D_LI_t/D_II_t*D_IL_t;
    D_LR=D_LR_t-D_LI_t/D_II_t*D_IR_t;
    D_RL=D_RL_t-D_RI_t/D_II_t*D_IL_t;
    D_RR=D_RR_t-D_RI_t/D_II_t*D_IR_t;
    
    % solving the polynomial eigenvalue problem:  
    [X,e]=polyeig(D_RL,(D_RR+D_LL),D_LR); 
    [Y,indices]=sort(real(log(e)));
    k_times_L(:,ind)=log(e(indices));
    phi(:,:,ind)=X(:,indices);
end

% Plotting results

scrsz = get(0,'ScreenSize');
width=550;
height=350;

figure('Position',[(scrsz(3)-width)/2 50 width height*2])
subplot(211)
hold on
plot3(omega/(2*pi),real(k_times_L),imag(k_times_L),'.','linewidth',2)
ylabel('Dimensionless real wavenumber k.L [-]')
% set(gca,'YTick',[-5*pi:pi:5*pi],'YGrid','on')
zlabel('Dimensionless imagenary wavenumber k.L [-]')
xlabel('Dimensionless frequency f/f_{natural beam} [-]')
set(gca,'ZTick',[-5*pi:pi/2:5*pi],'ZGrid','on')
title({['Euler-Bernouilli beam: E=',num2str(E/10^9), 'GPa, rho= ',num2str(rho),'kg/m^3','','','' ];...
    ['Mass ratio =',num2str(m_ratio),', f_{spring resonance} =',num2str(f_band_low_dimensionless),'.f_{natural beam}','Damping =',num2str(ksi),'.c_{critical}']})
grid on

subplot(212)
hold on
plot3(omega/(2*pi),real(k_times_L),imag(k_times_L),'.','linewidth',2)
ylabel('Dimensionless real wavenumber k.L [-]')
zlabel('Dimensionless imagenary wavenumber k.L [-]')
xlabel('Dimensionless frequency f/f_{natural beam} [-]')
set(gca,'ZTick',[-5*pi:pi/2:5*pi],'ZGrid','on')
grid on
view(0,0)

plot(omega/(2*pi),imag(k_times_L),'linewidth',2)

