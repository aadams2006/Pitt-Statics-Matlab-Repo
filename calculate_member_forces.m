% CALCULATE_MEMBER_FORCES.M
%   Use the 2D truss stiffness method to compute axial forces in each
%   bridge member. Update the node coordinates, connectivity, supports,
%   and applied loads in the input section to match your bridge model.

%% Input data ------------------------------------------------------------
nodeCoords = [ ...            % [x, y] coordinates (m)
    0.0, 0.0;  % Node 1 - left support
    4.0, 0.0;  % Node 2 - lower panel point
    8.0, 0.0;  % Node 3 - right support
    4.0, 3.0]; % Node 4 - upper panel point

E = 200e9;                     % Modulus of elasticity for steel (Pa)
A = 4.0e-4;                    % Cross-sectional area for all members (m^2)

members = struct( ...
    'name', {"BottomChord_L", "BottomChord_R", "Diagonal_L", "Diagonal_R", "TopChord"}, ...
    'nodes', {[1, 2], [2, 3], [1, 4], [4, 3], [4, 2]}, ...
    'area',  {A, A, A, A, A}, ...
    'E',     {E, E, E, E, E});

% Boundary conditions: restrained degrees of freedom (node, direction)
% direction: 1 = x, 2 = y
supports = [ ...
    1, 1;   % Node 1, ux = 0
    1, 2;   % Node 1, uy = 0
    3, 2];  % Node 3, uy = 0 (roller support)

% Applied loads [node, Fx (N), Fy (N)]
loads = [ ...
    2, 0.0, -30e3;   % Downward deck reaction at node 2
    3, 0.0, -30e3;   % Downward deck reaction at node 3
    4, 0.0, -10e3];  % Top panel point load

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
