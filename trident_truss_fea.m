%% FEA Project
%Authors ; Nathan FAGNOU
clc;
clear;
close all;

%% Data of the problem

% Variables
Nn = 4; % Number of nodes
Nr = 3; % Number of rods
L = 1; % Length of the rod
P = 1000; % Force applied on the system
alpha = pi/4; % Angle linked to P 

Module = ones(Nr, 1) * 200e9; % Young modulus E (Pa)
Section = ones(Nr, 1) * 1e-4;  % Section A (m^2)

Coord = [
    L, 0;  % Nod A 
    0, -L; % Nod B
    0, 0;  % Nod C
    0, L;  % Nod D
];

% Connectivity of elements
Connec = [
    1, 2; % Rod 1 (A-B)
    1, 3; % Rod 2 (A-C)
    1, 4; % Rod 3 (A-D)
];
% ---------------------------------------------------------------------

% P decomposition
P_X = P * cos(alpha);
P_Y = P * sin(alpha);

% Matrix orientation / length creation
Length = zeros(Nr, 1);
Orient = zeros(Nr, 1);

for i = 1:Nr
    n1 = Connec(i, 1);
    n2 = Connec(i, 2);

    x1 = Coord(n1, 1); y1 = Coord(n1, 2);
    x2 = Coord(n2, 1); y2 = Coord(n2, 2);

    dX = x2 - x1;
    dY = y2 - y1;

    Length(i) = sqrt(dX^2 + dY^2); % Using Pythagore
    Orient(i) = atan2(dY, dX);
end

%% Boundary conditions
% BC: Ux = 0 and Uy = 0 for nodes 2,3,4 (B, C, D)
BC_Ux = ones(Nn, 1);   % 1 means the node is fixed 
BC_Uy = ones(Nn, 1);   % same
VAL_Ux = zeros(Nn, 1); 
VAL_Uy = zeros(Nn, 1);

BC_Ux(1) = 0; BC_Uy(1) = 0;  % A is the only nod free  


% Force BCs
BC_Fx = zeros(Nn, 1);  % No external forces except for A (see after )
BC_Fy = zeros(Nn, 1);  % Same following Y 
VAL_Fx = zeros(Nn, 1);
VAL_Fy = zeros(Nn, 1);

% load P at node 1
BC_Fx(1) = 1; BC_Fy(1) = 1;
VAL_Fx(1) = P_X;
VAL_Fy(1) = P_Y;

%initialisation of vector F_global and q_global
N_dof_global = 2 * Nn;  % Definition of the number of degrees of freedom
Row_Column_elimin = zeros(N_dof_global,1);  % allow to set DOF for each nods
q_global = zeros(N_dof_global,1);  % Initialisation of displacement vector
F_global = zeros(N_dof_global,1);  % Same for the force 

for i = 1:Nn
    idx_ux = 2*i - 1; % 1 is X1, 3 is X3, etcc
    idx_uy = 2*i; % same for y
    Row_Column_elimin(idx_ux) = BC_Ux(i);  % applied to each nod
    Row_Column_elimin(idx_uy) = BC_Uy(i);
    q_global(idx_ux) = VAL_Ux(i);
    q_global(idx_uy) = VAL_Uy(i);
    F_global(idx_ux) = VAL_Fx(i);
    F_global(idx_uy) = VAL_Fy(i);
end

%% Definition of global stiffness matrix K 
N_dof_global = 2 * Nn; % Number of degrees of freedom
K = zeros(N_dof_global, N_dof_global);

for i = 1:Nr
    E = Module(i);
    A = Section(i);
    L_rod = Length(i);

    n1_idx = Connec(i, 1);
    n2_idx = Connec(i, 2);

    dX = Coord(n2_idx, 1) - Coord(n1_idx, 1);
    dY = Coord(n2_idx, 2) - Coord(n1_idx, 2);

    c = dX / L_rod;
    s = dY / L_rod;
    % Element stiffness matrix
    Ke = (E * A / L_rod) * [
        c^2,  c*s,  -c^2,  -c*s;
        c*s,  s^2,  -c*s,  -s^2;
       -c^2, -c*s,  c^2,   c*s;
       -c*s, -s^2,  c*s,   s^2;
    ];

    % Display element stiffness
    fprintf('\nElement %d: length = %.4f m, cos = %1.6f, sin = %1.6f\n', i, L_rod, c, s);
    fprintf('Ke (element %d):\n', i);
    disp(Ke);

    dof_n1 = [2*n1_idx - 1, 2*n1_idx];
    dof_n2 = [2*n2_idx - 1, 2*n2_idx];

    % Global DOF indices for this element
    dof_global = [dof_n1, dof_n2];
    for j = 1:4
        for k = 1:4
            K(dof_global(j), dof_global(k)) = K(dof_global(j), dof_global(k)) + Ke(j,k);
        end
    end
end

fprintf('\n Global stiffness matrix K \n');
disp(K);

%% Reduced stiffness matrix and reduced load vector  % suppression of fixed dof

dof_free = find(Row_Column_elimin == 0); % vector for nod of unknown displacment
dof_fixed = find(Row_Column_elimin == 1);
q_fixed = q_global(dof_fixed); 

N_BC_U = sum(Row_Column_elimin); % sum of fixed dof
N_dof_reduced = N_dof_global - N_BC_U; % sie of the system

if N_dof_reduced > 0
    K_red = K(dof_free, dof_free);
    fprintf('\n Reduced stiffness matrix K_red (kept DOFs only)\n');
    disp(K_red);

    F_free = F_global(dof_free);  % force vector for A
    K_fc = K(dof_free, dof_fixed);
    Correction_term = K_fc * q_fixed;
    F_red = F_free - Correction_term;
end
    fprintf('Reduced load vector F_red \n');
    disp(F_red');

%% Solve reduced system
if N_dof_reduced > 0
    q_red = K_red \ F_red; % to resolve our linear system
else
    q_red = [];
end

% Rebuild full displacement vector q
if N_dof_reduced > 0
    q_red = K_red \ F_red; % to resolve our linear system
else
    q_red = [];
end
% Rebuild full displacement vector q
q = zeros(N_dof_global, 1);
q(dof_free) = q_red; % assigned value for dof
q(dof_fixed) = q_fixed; % same thing for fixed dof
Nodal_Displacements = reshape(q, 2, Nn)'; % [Ux, Uy] per node

fprintf('\n---------------------------------------------------\n');
fprintf('Load and alpha\n');
fprintf('Load P = %.2f N, Angle alpha = %.2f deg\n', P, rad2deg(alpha));
fprintf('---------------------------------------------------\n\n');

fprintf('Nodal displacements (Ux, Uy) [m]:\n');
fprintf('-------------------------------------------------------------------\n');
for i = 1:Nn
    fprintf('Nœud %d | U_x%d = %1.10e (m) | U_y%d = %1.10e (m)\n', ...
            i, i, Nodal_Displacements(i,1), i, Nodal_Displacements(i,2));
end
fprintf('-------------------------------------------------------------------\n');

%% Reaction forces
R_global = K * q;
R_nodal = reshape(R_global, 2, Nn)';

fprintf('\nReaction forces at constrained DOFs and applied forces (Rx, Ry) [N]:\n');
fprintf('Node | Rx (N)            | Ry (N)\n');
fprintf('---------------------------------------------\n');
for i = 1:Nn
    if Row_Column_elimin(2*i-1) == 1 || Row_Column_elimin(2*i) == 1
        Reaction_X = R_nodal(i,1);
        Reaction_Y = R_nodal(i,2);

        if Row_Column_elimin(2*i-1) == 0
            Reaction_X = F_global(2*i-1);
        end
        if Row_Column_elimin(2*i) == 0
            Reaction_Y = F_global(2*i);
        end
        fprintf('%4d | %1.6f           | %1.6f\n', i, Reaction_X, Reaction_Y);
    end
end
fprintf('---------------------------------------------\n');

%% Visualization 
scale_factor = 1000;
Coord_deformed = Coord + scale_factor * Nodal_Displacements;
figure;
hold on;
grid on;
title(sprintf('Truss configuration (deformed scaled x %d)', scale_factor));
% Horizontal axis is Y
xlabel('Y coordinate (m)'); 
% Vertical axis is X
ylabel('X coordinate (m)'); 
% Use to invert axis 
set(gca, 'YDir', 'reverse');
% ------------------------------------------------------------------------------------------

% Plot initial configuration (y vs x)
for i = 1:Nr
    n1 = Connec(i,1);
    n2 = Connec(i,2);
    if i == 1
        h_initial = plot([Coord(n1,2), Coord(n2,2)], [Coord(n1,1), Coord(n2,1)], 'k--');
    else
        plot([Coord(n1,2), Coord(n2,2)], [Coord(n1,1), Coord(n2,1)], 'k--');
    end
end

% Plot deformed configuration
for i = 1:Nr
    n1 = Connec(i,1);
    n2 = Connec(i,2);
    if i == 1
        h_deformed = plot([Coord_deformed(n1,2), Coord_deformed(n2,2)], [Coord_deformed(n1,1), Coord_deformed(n2,1)], 'r-','LineWidth',1.5);
    else
        plot([Coord_deformed(n1,2), Coord_deformed(n2,2)], [Coord_deformed(n1,1), Coord(n2,1)], 'r-','LineWidth',1.5);
    end
end

plot(Coord(:,2), Coord(:,1), 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
text(Coord(:,2), Coord(:,1), arrayfun(@(k) sprintf(' %d',k), (1:Nn)','UniformOutput',false), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
plot(Coord_deformed(:,2), Coord_deformed(:,1), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
legend([h_initial, h_deformed], 'Initial', sprintf('Deformed (x%d)', scale_factor), 'Location', 'best');
axis equal;
hold off;