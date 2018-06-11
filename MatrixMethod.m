function [R] = MatrixMethod(L,theta,qy,qx,Ixx0,Iyy0,E,bcs)

%Author: Louis Kokee
%This function calculates the beam internal bending moment at each node.
%The calculations can be done for statically indeterminate straight beams
%and are done using the Stiffness Matrix method.
% clear all; close all; clc;


%% convert coordinate system
qz = qx;

% %%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%
% L = 1.0;                % m
% theta = 00;             % deg, rotation of coordinate frame from (x,y',z') to (x,y,z)
                      
Py = 0e3;               % N
Pz = 0e3;               % N  
% qy = 1.0e3;               % N/m
% qz = 1.0e3;               % N/m
                      
N = 100;              % Number of nodes, number will increase a bit during discretization

% %%%%%%%%%%%%%%%%% SECTION PROPERTIES %%%%%%%%%%%%%%%%%%%%%%
% % Moments of area
Iz = Ixx0;              % m^4
Iy = Iyy0;             % m^4

E = 1e5;            % Pa, Young's Modulus

%%%%%%%%%%%%%%%%% BOUNDARY CONDITIONS & FORCES %%%%%%%%%%%%%%%%%%%%%%
% For BC arrays:
% Legend ----> index0: x-location, index1: DoF, index2: value

% Each node has 4 DoF ordered as follows:
% Forces [Fy, Mz, Fz, My], Displacements [u_y, theta_z, u_z, theta_y]$
% Example: in the following arrays, if index2 = 2 it refers to Fz

% Applying a force at the tip of the structure
fext = [[L,1,Py];
        [L,3,Pz]];
        %[0.2,2,300]];

% Setting clamped BC at x = 0
bcs = [[0.0, 1, 0];
       [0.0, 2, 0];
       [0.0, 3, 0];
       [0.0, 4, 0];
       [L, 1, 0];
       [L, 2, 0];
       [L, 3, 0];
       [L, 4, 0]];

%%%%%%%%%%%%%%%%% PRE PROCESSING %%%%%%%%%%%%%%%%%%%%%%
% Put all critical x-locations in a sorted array with no duplicates
bc_loc = [bcs(:,1);fext(:,1)]; bc_loc = sort(unique(bc_loc(:,1)));
 
% Call discretization function
fprintf('\nDiscretizing...\n')
[xnodes, N] = discretize(L, bc_loc, N);

fprintf('New number of nodes:%i\n',N)

% Generate dictionary where key = x-location, value = corresponding node index
for i=1:length(bc_loc); idx(i)=find(xnodes==bc_loc(i)); end
id_xloc = containers.Map(bc_loc,idx);

% Keep track of fixed and free degree of freedoms
dof = 4;                             % DoFs
dofs_fixed = [];                     % array containing fixed DoFs
dofs_free = 1:dof*N;                 % array containing free DoFs

% Initialize matrices, 4 DoFs per node
fprintf('\nInitializing matrices...\n')
K = zeros(N*dof,N*dof); 
fvect = zeros(N*dof,1);
dvect = zeros(N*dof,1);

%% loads and boundary conditions
fprintf('\nApplying boundary conditions and external loads...\n')

% Applying discrete external loads to vector fvect
for ii = 1:length(fext(:,1))   
    fvect((id_xloc(fext(ii,1))-1)*dof + fext(ii,2)) = fext(ii,3);
end

% applying distributed external loads

%applying boundary conditions
for ii = 1:length(bcs(:,1))
	idx = (id_xloc(bcs(ii,1)) - 1)*dof + bcs(ii,2);
	dofs_fixed(ii) = idx;
% 	dofs_free(idx) = [];
	dvect(idx) = bcs(ii,3);
end
dofs_free(dofs_fixed) = [];


fprintf('\nLoading global stiffness matrix and force vector...\n')   
for iel =  1:N-1

	Le = xnodes(iel+1) - xnodes(iel);
	Kel = elStiffMatr(Le,E,Iz,Iy,theta);
    
    %load element stiffness matr in global stiffness matr
	K((iel-1)*dof + 1 : (iel+1)*dof,(iel-1)*dof + 1 : (iel+1)*dof) =...
        K((iel-1)*dof + 1 : (iel+1)*dof,(iel-1)*dof + 1 : (iel+1)*dof) + Kel;
        
    %effect of distributed force q
	fvect((iel-1)*dof + 1 : (iel+1)*dof) = ... 
        fvect((iel-1)*dof + 1 : (iel+1)*dof) + elNodeReaction(qy,qz,Le);
end
 
%% SOLVING 
fprintf('\nSolving system...\n')
K1 = K(dofs_free,dofs_free);
K12 = K(dofs_free,dofs_fixed);

fprintf('\nSolving for deflections...\n')
dvect(dofs_free) = K1 \ (fvect(dofs_free) - K12*dvect(dofs_fixed));

fprintf('\nSolving for reaction forces...\n')
R = bcs;            %array containing reaction forces, same structure as bcs
R(:,3) = (K(dofs_fixed,:)*dvect - fvect(dofs_fixed));


% %% convert coordinate system
% for i=length(R(:,1))
%     if R(i,2)==2
        


% figure
% subplot(2,1,1)
% plot(xnodes,1.0e3*dvect(1:4:end))
% xlabel('x [m]')
% ylabel('y deflection [mm]')
% 
% subplot(2,1,2)
% plot(xnodes,1.0e3*dvect(3:4:end))
% ylabel('z deflection [mm]')
% 
% % 3D Plot
% figure
% plot3(xnodes,1.0e3*dvect(1:4:end),1.0e3*dvect(3:4:end))
% grid on
% rotate3d on
% xlabel('x [m]')
% ylabel('y deflection [mm]')
% zlabel('z deflection [mm]')

end

