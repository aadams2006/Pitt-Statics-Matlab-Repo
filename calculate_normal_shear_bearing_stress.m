% CALCULATE_NORMAL_SHEAR_BEARING_STRESS.M
%   Compute normal, shear, and bearing stresses for bridge members.
%   Update the input data section with the forces and areas for the
%   connection or member you want to evaluate. Units are SI.

%% Input data ------------------------------------------------------------
% Positive axial force = tension, negative = compression (Newtons)
memberData = struct( ...
    'name',        {"TopChord", "BottomChord", "Diagonal"}, ...
    'axialForce',  {120e3,      -95e3,          60e3}, ...   % N
    'area',        {3.5e-4,     3.0e-4,         2.5e-4}, ... % m^2
    'shearForce',  {18e3,       12e3,           10e3}, ...   % N
    'bearingForce',{40e3,       35e3,           20e3}, ...   % N
    'bearingArea', {1.5e-3,     1.2e-3,         0.9e-3});    % m^2

%% Stress calculations ---------------------------------------------------
memberTable = struct2table(memberData);
memberTable.NormalStress_MPa  = normalStress([memberData.axialForce], [memberData.area]) * 1e-6;
memberTable.ShearStress_MPa   = shearStress([memberData.shearForce], [memberData.area]) * 1e-6;
memberTable.BearingStress_MPa = bearingStress([memberData.bearingForce], [memberData.bearingArea]) * 1e-6;

%% Output ----------------------------------------------------------------
disp("Stress summary (positive = tension / bearing / shear):");
disp(memberTable(:, ["name", "NormalStress_MPa", "ShearStress_MPa", "BearingStress_MPa"]));

%% Helper functions ------------------------------------------------------
function sigma = normalStress(P, A)
%NORMALSTRESS Compute the axial (normal) stress for an axial load.
    sigma = P ./ A;
end

function tau = shearStress(V, A)
%SHEARSTRESS Compute the average shear stress on a shear plane.
    tau = V ./ A;
end

function sigma_b = bearingStress(P, A_b)
%BEARINGSTRESS Compute the bearing stress at a connection interface.
    sigma_b = P ./ A_b;
end
