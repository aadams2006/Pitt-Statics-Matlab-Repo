% Simplified Fink Truss Analysis
% 200 lbf load at center, Poplar wood 9.525mm x 9.525mm
clear; clc;

%% UNIT CONVERSIONS
in_to_mm = 25.4;   % Convert inches to millimeters
lb_to_N = 4.448;   % Convert pounds-force to Newtons

%% GEOMETRY
% Joint coordinates [x, y] in inches (converted to mm below)
% These define the 9 connection points of the truss after adding two new joints
base_joints_in = [
    -8.5,  0.0;   % 1 - left support (pin - fixed in x and y)
    -2.93, 0.0;   % 2 - bottom chord connection point
     2.93, 0.0;   % 3 - bottom chord connection point
     8.5,  0.0;   % 4 - right support (roller - fixed in y only)
    -4.4,  2.54;  % 5 - left upper joint
     0.0,  5.08;  % 6 - apex (top center)
     4.4,  2.54;  % 7 - right upper joint
];

% Determine optimal (45°) brace connections for new members
target_angle_deg = 45;  % Symmetric 45° brace to minimize member length and provide efficient load path
target_angle_rad = deg2rad(target_angle_deg);

% New joint on member 1-2 connected to joint 5
dy_left_in = base_joints_in(5, 2);                % Vertical offset from joint 5 to bottom chord (in)
dx_left_in = dy_left_in / tan(target_angle_rad);
joint8_in = [base_joints_in(5, 1) - dx_left_in, 0];  % Joint 8 - brace connection on member 1-2

% New joint on member 3-4 connected to joint 7
dy_right_in = base_joints_in(7, 2);
dx_right_in = dy_right_in / tan(target_angle_rad);
joint9_in = [base_joints_in(7, 1) + dx_right_in, 0];  % Joint 9 - brace connection on member 3-4

% Register joint indices before appending (clarifies the new joint numbering)
joint8_index = size(base_joints_in, 1) + 1;
joint9_index = joint8_index + 1;

% Append new joints and convert to millimeters
joints_in = [
    base_joints_in;
    joint8_in;
    joint9_in;
];
joints = joints_in * in_to_mm;

fprintf('Joint %d placed on member 1-2 at (%.3f, %.3f) in\n', joint8_index, joint8_in(1), joint8_in(2));
fprintf('Joint %d placed on member 3-4 at (%.3f, %.3f) in\n', joint9_index, joint9_in(1), joint9_in(2));

% Report the resulting brace angles (acute angle relative to the bottom chord)
angle_5_to_8 = abs(atand((joints(5,2) - joints(joint8_index,2)) / ...
                        (joints(joint8_index,1) - joints(5,1))));
angle_7_to_9 = abs(atand((joints(7,2) - joints(joint9_index,2)) / ...
                        (joints(joint9_index,1) - joints(7,1))));
fprintf('New member 5-%d angle to horizontal: %.2f degrees\n', joint8_index, angle_5_to_8);
fprintf('New member 7-%d angle to horizontal: %.2f degrees\n', joint9_index, angle_7_to_9);

% Members [start, end]
% Note: Bottom chord (members 1-5) is physically one continuous beam
% but modeled as 5 segments for truss analysis (including new brace joints)
members = [
    1, joint8_index;    % 1  - Bottom chord (left outer segment)
    joint8_index, 2;    % 2  - Bottom chord (left inner segment)
    2, 3;               % 3  - Bottom chord (center segment)
    3, joint9_index;    % 4  - Bottom chord (right inner segment)
    joint9_index, 4;    % 5  - Bottom chord (right outer segment)
    1, 5;               % 6  - Left diagonal
    2, 5;               % 7  - Left vertical
    5, 6;               % 8  - Left diagonal
    2, 6;               % 9  - Left diagonal
    3, 6;               % 10 - Right diagonal
    6, 7;               % 11 - Right diagonal
    3, 7;               % 12 - Right vertical
    4, 7;               % 13 - Right diagonal
    5, joint8_index;    % 14 - New left brace
    7, joint9_index;    % 15 - New right brace
];

%% MATERIAL PROPERTIES
width = 9.525;                            % mm (3/8")
height = 9.525;                           % mm (3/8")
A = width * height;                       % Cross-sectional area = 90.7 mm²
E = 8400;                                 % Young's modulus = 8400 MPa (compression)

% Allowable stresses for Poplar wood (MPa)
sigma_tension_allow = 75;                 % Maximum tension stress
sigma_compression_allow = 38;             % Maximum compression stress
tau_shear_allow = 45;                     % Maximum shear stress
sigma_bearing_allow = 38;                 % Maximum bearing stress (crushing)

% Bolt/Pin connection properties
bolt_diameter = 4.166;                    % mm (0.164")
A_bolt = pi * (bolt_diameter/2)^2;        % Bolt cross-sectional area (mm²)
A_bearing = bolt_diameter * width;        % Bearing area = bolt diameter × wood thickness

%% LOADING & SUPPORTS
% 250 lbf total load split equally between nodes 2 and 3 (center of bridge)
loads = [2, 0, -100; 3, 0, -100];  % [joint#, Fx, Fy] in lbf
% Fixed DOFs: DOF 1 = joint1_x, DOF 2 = joint1_y, DOF 8 = joint4_y
% This creates: pin support at joint 1 (can't move) and roller at joint 4 (can slide horizontally)
fixed_dofs = [1, 2, 8];

%% SOLVE FOR DISPLACEMENTS
n_joints = size(joints, 1);    % Number of joints = 9
n_members = size(members, 1);  % Number of members = 15
n_dof = 2 * n_joints;          % Total degrees of freedom = 18 (2 per joint: x and y)

% Initialize global stiffness matrix (18×18) and force vector (18×1)
K_global = zeros(n_dof);
F = zeros(n_dof, 1);

% Assemble global stiffness matrix
% This loop builds the structural stiffness matrix by adding each member's contribution
for i = 1:n_members
    j1 = members(i, 1);  % Start joint number
    j2 = members(i, 2);  % End joint number
    
    % Calculate member geometry
    dx = joints(j2,1) - joints(j1,1);  % Change in x
    dy = joints(j2,2) - joints(j1,2);  % Change in y
    L = sqrt(dx^2 + dy^2);             % Member length
    c = dx/L;  % Cosine of angle (direction cosine in x)
    s = dy/L;  % Sine of angle (direction cosine in y)
    
    % Member axial stiffness: k = (Area × Young's Modulus) / Length
    k = (A * E) / L;
    
    % Transformation matrix: converts local (along member) to global (x,y) coordinates
    T = [c, s, 0, 0;
         0, 0, c, s];
    
    % Local stiffness matrix (tension-compression only)
    k_local = k * [1, -1; -1, 1];
    
    % Transform to global coordinates
    k_global = T' * k_local * T;
    
    % Add member stiffness to global matrix at appropriate DOFs
    % DOFs for joints: [j1_x, j1_y, j2_x, j2_y]
    dofs = [2*j1-1, 2*j1, 2*j2-1, 2*j2];
    K_global(dofs, dofs) = K_global(dofs, dofs) + k_global;
end

% Apply loads to force vector
for i = 1:size(loads, 1)
    joint = loads(i, 1);
    F(2*joint-1) = F(2*joint-1) + loads(i, 2) * lb_to_N;  % Fx (convert lbf to N)
    F(2*joint) = F(2*joint) + loads(i, 3) * lb_to_N;      % Fy (convert lbf to N)
end

% Solve for displacements
% Remove rows/columns for fixed DOFs (supports), then solve K*U = F
free_dofs = setdiff(1:n_dof, fixed_dofs);  % Find which DOFs can move
U = zeros(n_dof, 1);                        % Initialize displacement vector
U(free_dofs) = K_global(free_dofs, free_dofs) \ F(free_dofs);  % Solve for free DOFs

%% ========================================
%% 1. MEMBER FORCES & NORMAL STRESSES
%% ========================================
fprintf('\n===============================================================\n');
fprintf('1. MEMBER FORCES & NORMAL STRESSES\n');
fprintf('===============================================================\n');
fprintf('Member  Force(lbf) Type  Normal Stress(MPa)  Allowable  Status\n');
fprintf('------  ---------- ----  ------------------  ---------  ------\n');

member_forces = zeros(n_members, 1);
for i = 1:n_members
    j1 = members(i, 1);
    j2 = members(i, 2);
    
    % Member geometry (same as before)
    dx = joints(j2,1) - joints(j1,1);
    dy = joints(j2,2) - joints(j1,2);
    L = sqrt(dx^2 + dy^2);
    c = dx/L;
    s = dy/L;
    
    % Get displacements at each end of member
    u1 = U(2*j1-1); v1 = U(2*j1);  % Joint 1 displacements (x, y)
    u2 = U(2*j2-1); v2 = U(2*j2);  % Joint 2 displacements (x, y)
    
    % Calculate axial force: F = (E*A/L) × elongation
    % Elongation = projection of displacement difference onto member axis
    force_N = (E * A / L) * (c*(u2-u1) + s*(v2-v1));
    force_lbf = force_N / lb_to_N;  % Convert to lbf for display
    
    % Normal stress: σ = Force / Area
    stress = force_N / A;
    
    member_forces(i) = force_lbf;
    
    % Determine type and check against appropriate allowable stress
    if force_N > 0.01
        type = 'T';  % Tension (positive force = pulling)
        allowable = sigma_tension_allow;
    elseif force_N < -0.01
        type = 'C';  % Compression (negative force = pushing)
        allowable = sigma_compression_allow;
        stress = abs(stress);  % Use absolute value for comparison
    else
        type = 'Z';  % Zero force (very small)
        allowable = sigma_compression_allow;
        stress = 0;
    end
    
    % Check if stress is within allowable limit
    status = iif(stress <= allowable, 'PASS', 'FAIL');
    
    fprintf('%3d     %10.2f %s     %13.2f       %9.2f  %s\n', ...
            i, force_lbf, type, stress, allowable, status);
end
fprintf('\nNote: T=Tension, C=Compression, Z=Zero\n');

%% ========================================
%% 2. SHEAR STRESS AT JOINTS (Double Shear)
%% ========================================
fprintf('\n===============================================================\n');
fprintf('2. SHEAR STRESS AT JOINTS (Double Shear)\n');
fprintf('===============================================================\n');
fprintf('Joint  Max Force(lbf)  Shear Stress(MPa)  Allowable  Status\n');
fprintf('-----  -------------  -----------------  ---------  ------\n');

% Find maximum force at each joint from all connected members
joint_forces = zeros(n_joints, 1);
for i = 1:n_members
    j1 = members(i, 1);
    j2 = members(i, 2);
    force_mag = abs(member_forces(i));
    % Store the largest force acting at each joint
    joint_forces(j1) = max(joint_forces(j1), force_mag);
    joint_forces(j2) = max(joint_forces(j2), force_mag);
end

for i = 1:n_joints
    % Double shear: bolt is cut at TWO cross-sections
    % τ = Force / (2 × bolt area)
    % Factor of 2 because load is shared between two shear planes
    tau_shear = (joint_forces(i) * lb_to_N) / (2 * A_bolt);
    
    status = iif(tau_shear <= tau_shear_allow, 'PASS', 'FAIL');
    fprintf('%3d    %13.2f  %17.2f  %9.2f  %s\n', ...
            i, joint_forces(i), tau_shear, tau_shear_allow, status);
end

%% ========================================
%% 3. BEARING STRESS AT JOINTS
%% ========================================
fprintf('\n===============================================================\n');
fprintf('3. BEARING STRESS AT JOINTS\n');
fprintf('===============================================================\n');
fprintf('Joint  Max Force(lbf)  Bearing Stress(MPa)  Allowable  Status\n');
fprintf('-----  -------------  -------------------  ---------  ------\n');

for i = 1:n_joints
    % Bearing stress: σ_bearing = Force / (bolt diameter × wood thickness)
    % This is the crushing stress on the wood where the bolt contacts it
    sigma_bearing = (joint_forces(i) * lb_to_N) / A_bearing;
    
    status = iif(sigma_bearing <= sigma_bearing_allow, 'PASS', 'FAIL');
    fprintf('%3d    %13.2f  %19.2f  %9.2f  %s\n', ...
            i, joint_forces(i), sigma_bearing, sigma_bearing_allow, status);
end

%% ========================================
%% 4. BOTTOM CHORD BEAM SHEAR STRESS
%% ========================================
fprintf('\n===============================================================\n');
fprintf('4. BOTTOM CHORD BEAM SHEAR STRESS (Continuous Bottom Beam)\n');
fprintf('===============================================================\n');

% Bottom chord acts as a continuous beam from joint 1 to joint 4
beam_span = abs(joints(4,1) - joints(1,1));  % Total span (mm)
total_load_N = 250 * lb_to_N;                 % Total load in N

% For a simply supported beam with center point load P:
% Maximum shear occurs at supports and equals P/2
V_max = total_load_N / 2;  % Max shear force (N)

% Shear stress for rectangular cross-section: τ = 1.5 × V / Area
% Factor of 1.5 accounts for parabolic shear distribution
tau_beam = 1.5 * V_max / A;

fprintf('Bottom chord span:      %.2f mm\n', beam_span);
fprintf('Total load:             %.2f lbf\n', 250);
fprintf('Max shear force:        %.2f lbf (at supports)\n', V_max / lb_to_N);
fprintf('Max shear stress:       %.2f MPa\n', tau_beam);
fprintf('Allowable shear:        %.2f MPa\n', tau_shear_allow);
fprintf('Status:                 %s\n', iif(tau_beam <= tau_shear_allow, 'PASS', 'FAIL'));

%% SUMMARY
fprintf('\n===============================================================\n');
fprintf('SUMMARY - ALL CHECKS\n');
fprintf('===============================================================\n');

% Check if all members pass their respective stress checks
all_normal_pass = all(abs(member_forces(member_forces>0)*lb_to_N/A) <= sigma_tension_allow) && ...
                  all(abs(member_forces(member_forces<0)*lb_to_N/A) <= sigma_compression_allow);
all_shear_pass = all((joint_forces*lb_to_N./(2*A_bolt)) <= tau_shear_allow);
all_bearing_pass = all((joint_forces*lb_to_N./A_bearing) <= sigma_bearing_allow);
cross_beam_pass = tau_beam <= tau_shear_allow;

fprintf('Normal stresses:        %s\n', iif(all_normal_pass, 'PASS', 'FAIL'));
fprintf('Shear at joints:        %s\n', iif(all_shear_pass, 'PASS', 'FAIL'));
fprintf('Bearing at joints:      %s\n', iif(all_bearing_pass, 'PASS', 'FAIL'));
fprintf('Bottom chord shear:     %s\n', iif(cross_beam_pass, 'PASS', 'FAIL'));
fprintf('\nOVERALL:                %s\n', ...
    iif(all_normal_pass && all_shear_pass && all_bearing_pass && cross_beam_pass, ...
        'ALL CHECKS PASS', 'DESIGN NEEDS REVIEW'));
fprintf('===============================================================\n');

%% Helper function
% Simple if-then-else function for returning strings
function result = iif(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end