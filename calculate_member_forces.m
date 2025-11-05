% CALCULATE_MEMBER_FORCES.M
%   Use the 2D truss stiffness method to compute axial forces in each
%   bridge member. Update the node coordinates, connectivity, supports,
%   and applied loads in the input section to match the bridge model.
%   This script is configured for a Fink Truss.

%% Input data ------------------------------------------------------------
nodeCoords = [ ...            % [x, y] coordinates (m)
    0.0, 0.0;  % Node 1 - left support
    5.87, 0.0;  % Node 2
    8.0, 0.0;  % Node 3
    11.73, 0.0; % Node 4
    16.0, 0.0; % Node 5 - right support
    4.0, 2.0;  % Node 6
    8.0, 4.0;  % Node 7
    12.0, 2.0; % Node 8
    8.0, 2.0]; % Node 9

E = 200e9;                     % Modulus of elasticity for steel (Pa)
A = 4.0e-4;                    % Cross-sectional area for all members (m^2)

members = struct( ...
    'name', {"BottomChord1", "BottomChord2", "BottomChord3", "BottomChord4", "TopChord1", "TopChord2", "TopChord3", "TopChord4", "Web1", "Web2", "Web3", "Web4", "Web5", "Web6", "Web7", "Web8"}, ...
    'nodes', {[1, 2], [2, 3], [3, 4], [4, 5], [1, 6], [6, 7], [7, 8], [8, 5], [2, 6], [2, 9], [9, 6], [9, 7], [3, 9], [4, 8], [4, 9], [9,8]}, ...
    'area',  {A, A, A, A, A, A, A, A, A, A, A, A, A, A, A, A}, ...
    'E',     {E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E});

% Boundary conditions: restrained degrees of freedom (node, direction)
% direction: 1 = x, 2 = y
supports = [ ...
    1, 1;   % Node 1, ux = 0
    1, 2;   % Node 1, uy = 0
    5, 2];  % Node 5, uy = 0 (roller support)

% Applied loads [node, Fx (N), Fy (N)]
loads = [ ...
    2, 0.0, -10e3;
    3, 0.0, -10e3;
    4, 0.0, -10e3;
    6, 0.0, -5e3;
    7, 0.0, -5e3;
    8, 0.0, -5e3];

%% Assemble global stiffness matrix -------------------------------------
numNodes = size(nodeCoords, 1);
numDof = numNodes * 2;
K = zeros(numDof);

for m = 1:numel(members)
    n1 = members(m).nodes(1);
    n2 = members(m).nodes(2);
    x1 = nodeCoords(n1, 1); y1 = nodeCoords(n1, 2);
    x2 = nodeCoords(n2, 1); y2 = nodeCoords(n2, 2);

    L = hypot(x2 - x1, y2 - y1);
    c = (x2 - x1) / L;
    s = (y2 - y1) / L;

    kLocal = members(m).E * members(m).area / L * ...
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
forces = zeros(numel(members), 1);
for m = 1:numel(members)
    n1 = members(m).nodes(1);
    n2 = members(m).nodes(2);
    x1 = nodeCoords(n1, 1); y1 = nodeCoords(n1, 2);
    x2 = nodeCoords(n2, 1); y2 = nodeCoords(n2, 2);

    L = hypot(x2 - x1, y2 - y1);
    c = (x2 - x1) / L;
    s = (y2 - y1) / L;

    dofIdx = [dofIndex(n1, 1), dofIndex(n1, 2), dofIndex(n2, 1), dofIndex(n2, 2)];
    uElem = U(dofIdx);
    axialStrain = [-c, -s, c, s] * uElem / L;
    forces(m) = members(m).E * members(m).area * axialStrain; % Positive = tension
end

memberTable = struct2table(members);
memberTable.Force_kN = forces * 1e-3;

%% Output ----------------------------------------------------------------
disp("Member axial forces (positive = tension):");
disp(memberTable(:, ["name", "Force_kN"]));

disp("Node displacements (mm):");
disp(reshape(U, 2, numNodes)' * 1e3);

%% Helper function -------------------------------------------------------
function idx = dofIndex(node, direction)
%DOFINDEX Map a (node, direction) pair to a global degree-of-freedom index.
    idx = 2 * (node - 1) + direction;
end
