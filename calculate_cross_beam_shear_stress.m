% CALCULATE_CROSS_BEAM_SHEAR_STRESS.M
%   Estimate the shear stress distribution in a cross beam using the
%   elastic beam shear formula (VQ/Ib). Update the section geometry and
%   applied shear force to match your bridge design.

%% Input data ------------------------------------------------------------
V = 85e3;  % Shear force acting on the cross beam (N)

section = struct();
section.flangeWidth     = 0.18;  % Width of flanges (m)
section.flangeThickness = 0.018; % Thickness of each flange (m)
section.webThickness    = 0.010; % Web thickness (m)
section.webHeight       = 0.380; % Clear web height between flanges (m)

%% Section properties ----------------------------------------------------
H = 2 * section.flangeThickness + section.webHeight; % Overall depth (m)
I = (section.flangeWidth * H^3 - ...
    (section.flangeWidth - section.webThickness) * section.webHeight^3) / 12; % m^4

A_flange = section.flangeWidth * section.flangeThickness;
A_web = section.webThickness * section.webHeight;

% First moment of area about the neutral axis for the top p
ortion of the
% section (used to find maximum shear stress in the web at the neutral axis).
Q_na = A_flange * (section.webHeight / 2 + section.flangeThickness / 2) + ...
       (0.5 * A_web) * (section.webHeight / 4);

% Shear stress at the neutral axis (maximum in the web)
tau_max = V * Q_na / (I * section.webThickness);

% Average shear stress in the web for comparison
area_web = section.webThickness * section.webHeight;
tau_avg = V / area_web;

%% Output ----------------------------------------------------------------
fprintf('Cross beam shear analysis\n');
fprintf('-------------------------\n');
fprintf('Applied shear force:        %8.1f kN\n', V * 1e-3);
fprintf('Section moment of inertia:  %8.4e m^4\n', I);
fprintf('Max shear stress (web NA):  %8.2f MPa\n', tau_max * 1e-6);
fprintf('Average web shear stress:   %8.2f MPa\n', tau_avg * 1e-6);

%% Notes -----------------------------------------------------------------
% tau_max is computed at the neutral axis using the elastic beam theory.
% Modify Q_na if you need the shear stress at another location within the
% cross-section (e.g., flange/web interface).
