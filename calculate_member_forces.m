% CALCULATE_MEMBER_FORCES.M
%   Use the 2D truss stiffness method to compute axial forces in each
%   bridge member. Update the node coordinates, connectivity, supports,
%   and applied loads in the input section to match the bridge model.
%   This script is configured for a Fink Truss.

%% Input data ------------------------------------------------------------
% Node 3 is located at the origin on the bottom chord to match the bridge
% coordinate system used in the design models.
nodeCoords = [ ...            % [x, y] coordinates (m)
    -6.0, 0.0;   % Node 1  - left support
    -2.0, 0.0;   % Node 2
     0.0, 0.0;   % Node 3  - bottom chord center (origin)
     6.0, 0.0;   % Node 4  - right support
    -3.0, 3.0;   % Node 5
     0.0, 5.0;   % Node 6 - apex
     3.0, 3.0];  % Node 7

E = 200e9;                     % Modulus of elasticity for steel (Pa)
A = 4.0e-4;                    % Cross-sectional area for all members (m^2)

memberNames = {
    "BottomChord1", "BottomChord2", "BottomChord3", ...
    "TopChord1", "TopChord2", ...
    "LeftEndPost", "RightEndPost", ...
    "LeftDiagonalLower", "LeftDiagonalUpper", ...
    "RightDiagonalUpper", "RightDiagonalLower"};

memberConnectivity = [ ...
    1, 2;  % Bottom chord
    2, 3;
    3, 4;
    5, 6;  % Top chord
    6, 7;
    1, 5;  % End posts
    4, 7;
    2, 5;  % Left lower diagonal
    2, 6;  % Left upper diagonal
    3, 6;  % Right upper diagonal
    3, 7]; % Right lower diagonal

numMembers = size(memberConnectivity, 1);
memberArea = A * ones(numMembers, 1);
memberE = E * ones(numMembers, 1);

% Boundary conditions: restrained degrees of freedom (node, direction)
% direction: 1 = x, 2 = y
supports = [ ...
    1, 1;   % Node 1, ux = 0
    1, 2;   % Node 1, uy = 0
    4, 2];  % Node 4, uy = 0 (roller support)

% Applied loads [node, Fx (N), Fy (N)]
loads = [ ...
    2, 0.0, -12e3;
    3, 0.0, -12e3];

%% Assemble global stiffness matrix -------------------------------------
numNodes = size(nodeCoords, 1);
numDof = numNodes * 2;
K = zeros(numDof);

for m = 1:numMembers
    n1 = memberConnectivity(m, 1);
    n2 = memberConnectivity(m, 2);
    x1 = nodeCoords(n1, 1); y1 = nodeCoords(n1, 2);
    x2 = nodeCoords(n2, 1); y2 = nodeCoords(n2, 2);

    L = hypot(x2 - x1, y2 - y1);
    c = (x2 - x1) / L;
    s = (y2 - y1) / L;

    kLocal = memberE(m) * memberArea(m) / L * ...
        [ c^2,  c*s, -c^2, -c*s; ...
          c*s,  s^2, -c*s, -s^2; ...
         -c^2, -c*s,  c^2,  c*s; ...
         -c*s, -s^2,  c*s,  s^2];

    dofIdx = [dofIndex(n1, 1), dofIndex(n1, 2), dofIndex(n2, 1), dofIndex(n2, 2)];
    K(dofIdx, dofIdx) = K(dofIdx, dofIdx) + kLocal;
end

%% Apply loads -----------------------------------------------------------
F = zeros(numDof, 1);
for i = 1:size(loads, 1)
    node = loads(i, 1);
    F(dofIndex(node, 1)) = F(dofIndex(node, 1)) + loads(i, 2);
    F(dofIndex(node, 2)) = F(dofIndex(node, 2)) + loads(i, 3);
end

%% Apply boundary conditions --------------------------------------------
constrainedDof = unique(dofIndex(supports(:, 1), supports(:, 2)));
freeDof = setdiff(1:numDof, constrainedDof);

U = zeros(numDof, 1);
U(freeDof) = K(freeDof, freeDof) \ F(freeDof);

%% Member force recovery -------------------------------------------------
forces = zeros(numMembers, 1);
memberLength = zeros(numMembers, 1);
for m = 1:numMembers
    n1 = memberConnectivity(m, 1);
    n2 = memberConnectivity(m, 2);
    x1 = nodeCoords(n1, 1); y1 = nodeCoords(n1, 2);
    x2 = nodeCoords(n2, 1); y2 = nodeCoords(n2, 2);

    L = hypot(x2 - x1, y2 - y1);
    memberLength(m) = L;
    c = (x2 - x1) / L;
    s = (y2 - y1) / L;

    dofIdx = [dofIndex(n1, 1), dofIndex(n1, 2), dofIndex(n2, 1), dofIndex(n2, 2)];
    uElem = U(dofIdx);
    axialStrain = [-c, -s, c, s] * uElem / L;
    forces(m) = memberE(m) * memberArea(m) * axialStrain; % Positive = tension
end

memberTable = table( ...
    memberNames', num2cell(memberConnectivity, 2), memberLength, forces * 1e-3, ...
    forces ./ memberArea * 1e-6, ...
    'VariableNames', {"Name", "Nodes", "Length_m", "Force_kN", "Stress_MPa"});

%% Output ----------------------------------------------------------------
disp("Member axial forces and stresses (positive = tension):");
disp(memberTable);

nodeDisplacements = reshape(U, 2, numNodes)' * 1e3;
nodeTable = table((1:numNodes)', nodeCoords(:, 1), nodeCoords(:, 2), ...
    nodeDisplacements(:, 1), nodeDisplacements(:, 2), ...
    'VariableNames', {"Node", "X_m", "Y_m", "Ux_mm", "Uy_mm"});

disp("Node displacements (mm):");
disp(nodeTable);

%% Helper function -------------------------------------------------------
function idx = dofIndex(node, direction)
%DOFINDEX Map a (node, direction) pair to a global degree-of-freedom index.
    idx = 2 * (node - 1) + direction;
end
