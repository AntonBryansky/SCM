clearvars;
close all;
%% Layers properties
% Properties of the base material (bi-directional CFRP T300 woven)
% p - density, kg/m^3
% E_s - Young modulus, Pa
% G_s - Shear modulus, Pa
% nu_s - Poisson's ratio
% b_s - Laminate thickness, m
p_s = 1437;
E_1 = 54.50*10^9;
E_2 = E_1;
E_3 = 8.4*10^9;
G_12 = 3.07*10^9;
G_13 = 4.2*10^9;
G_23 = G_13;
nu_12 = 0.0638;
nu_13 = 0.25;
nu_23 = nu_13;
b_s = 1e-3;

% Characteristic sizes of the honeycomb core
% t - thickness of the base material, m
% l - cell wall length, m
% h - cell wall length, m
% b - cell wall height, m
% theta - angle of the cell wall, pi
t = 0.21*10^-3;
l = 5*10^-3;
h = 5*10^-3;
b = 7e-3;
theta = 1/6;

% Calculate the homogenized elastic properties of the honeycomb core
[p_c, C_c] = SM_HC(t, l, h, b, theta, p_s, E_1, G_12, nu_12);
% Calculate the elastic properties of the skin sheet
C_s = SM_orthotropic(E_1, E_2, E_3, G_12, G_13, G_23, nu_12, nu_13, nu_23);

%% General parameters
% Mechanical properties of layers
densities = [p_s p_c p_s]; % density, kg/m^3
Cs = {C_s, C_c, C_s}; % stiffness matrices, Pa
layerOrientations = [0 0 0]; % orientation angles of the layers, degree
layerThicknesses = [b_s b b_s]; % thicknesses of the layers, m
% Differential matrices parameters, number of collocation points
N1 = 21;
N2 = 21;
N3 = 21;
Ns = [N1 N2 N3];
% This value is used instead of zero to avoid computational errors
z = 1e-15; % Error of calculation, depends on the number of collocation points

numOfLayers = length(layerThicknesses); % number of layers
propagationAngle = 0; % Wave propagation direction relative to layer orientation angles, degree

FreqLimit = 300;
WNLimit = 8000;
WNAmount = 5000;
wavenumber = 0:(WNLimit/WNAmount):WNLimit; wavenumber(1) = 1e-6;

is_orth = isOrthotropic(layerOrientations, propagationAngle); % check if the material ortothropic
is_sym = isSymmetric(densities, Cs, layerOrientations, layerThicknesses); % check if the layup is symmetric

% These options enable plotting phase and group velocities plotts, respectivelly
plot_GV = 1; % Group velocities plots
plot_PV = 1; % Phase velocities plots
% Additional options for Phase velocities plotting
plot_PV_all_modes = 1;
plot_PV_all_modes_separated = 1;

plot_PV_separated_A_S_modes = 1;
plot_PV_separated_L_SH_modes = 1;
plot_PV_L_separated_modes = 1;
plot_PV_SH_separated_modes = 1;

% This option enables the boundary condition collocations scheme:
% 1 - Quintanilla's scheme
% 2 - Modifiend Quintanilla's scheme
% 3 - Original Mekkaoui's scheme (doesn't work)
% 4 - Modified Mekkaoui's scheme
BC_scheme = 2;

% This option enables balancing algorhythm
% 0 - No balancing (numerically instable for multilayer problem)
% 1 - Matlab built-in balancing (only the eigenvalues)
% 2 - Ward algorithm (both eigenvectors and eigenvalues)
Balance_scheme = 2;

%% Transform stiffness matrices relatively to the propagation angle
beta = layerOrientations - propagationAngle; % angle between the main direction and propagation angly, usually equals to zero, degree
c = cell(1, numOfLayers); % Set transformed matrix
for i = 1:numOfLayers
    c{1, i} = transformStiffnessMatrix(Cs{i}, beta(i));
end

%% Calculation
% Independent constants
% Chebyshev differentiation matrices for each layer
DM = cell(1, numOfLayers);
DM2 = cell(1, numOfLayers);
Xs = cell(1, numOfLayers);
Is = cell(1, numOfLayers);
for i = 1:numOfLayers
    [x, D] = chebdif(Ns(i), 1);
    X_min = (sum(layerThicknesses)/2) - sum(layerThicknesses(1:(i)));
    X_max = (sum(layerThicknesses)/2) - sum(layerThicknesses(1:(i-1)));
    X_c = 0.5*(X_min + X_max);
    X = 0.5*(X_max - X_min)*x + X_c;
    Xs{1, i} = X;
    buff = D(:,:,1).*(2/(X_max - X_min));
    DM{1, i} = buff;
    DM2{1, i} = buff^2;
    Is{1, i} = eye(Ns(i));
end

% M matrices are independent of wavenumber
M = cell(1, numOfLayers);
for i = 1:numOfLayers
    M{1,i} = getMmatrix(Ns(i), densities(i));
end

% Output phase velocities
velocities = NaN(length(wavenumber),3*sum(Ns));
freqs = NaN(length(wavenumber), 3*sum(Ns));

Lwave = NaN(length(wavenumber), 3*sum(Ns));
SHwave = NaN(length(wavenumber), 3*sum(Ns));

Sym = NaN(length(wavenumber), 3*sum(Ns));
Asym = NaN(length(wavenumber), 3*sum(Ns));

SM = NaN(length(wavenumber), 3*sum(Ns));
ASM = NaN(length(wavenumber), 3*sum(Ns));
SHSM = NaN(length(wavenumber), 3*sum(Ns));
SHASM = NaN(length(wavenumber), 3*sum(Ns));

% Calculation dispersion curves
f = waitbar(0, 'Starting calculations...');
timestart = tic;
for i = 1:length(wavenumber)
    progr = i/length(wavenumber);
    waitbar(progr, f, strcat(num2str(round(progr*100, 3)), '%'))
    % calculating L, M and S matrices
    k = wavenumber(i);
    L = cell(1, numOfLayers);
    S = cell(1, numOfLayers);
    for n = 1:numOfLayers
        L{1,n} = getLmatrix(c{n}, k, Is{1, n}, DM{1, n}, DM2{1, n});
        S{1,n} = getSmatrix(c{n}, k, Is{1, n}, DM{1, n});
    end
    
    L_BC = zeros(3*sum(Ns), 3*sum(Ns));
    M_BC = zeros(3*sum(Ns), 3*sum(Ns));

    % M matrix boundary conditions
    for n = 1:numOfLayers
        M_BC((3*sum(Ns(1:(n-1)))+1):3*sum(Ns(1:n)), (3*sum(Ns(1:(n-1)))+1):3*sum(Ns(1:n))) = getMBC(M{1, n}, Ns(n));
    end
    
    % L matrix boundary conditions
    switch BC_scheme
        case 1 % Quintanilla's scheme of boundary conditions
            for n = 1:numOfLayers
                if n == 1 % First layer
                    NL = Ns(n); NR = Ns(n + 1);
                    buff_L = L{1, n};
                    buff_R = zeros(3*NL, 3*NR);
                    % External BC
                    buff_L(1, :) = S{1, n}(1, :);
                    buff_L((NL), :) = S{1, n}((NL+1), :);
                    buff_L((NL+1), :) = S{1, n}((2*NL+1), :);
                    % Continuity BC
                    % 1st layer
                    buff_L((2*NL), :) = S{1, n}(NL, :);
                    buff_L((2*NL+1), :) = S{1, n}((2*NL), :);
                    buff_L((3*NL), :) = S{1, n}((3*NL), :);
                    % 2nd layer
                    buff_R((2*NL), :) = -S{1, (n+1)}((1), :);
                    buff_R((2*NL+1), :) = -S{1, (n+1)}((NR+1), :);
                    buff_R((3*NL), :) = -S{1, (n+1)}((2*NR+1), :);
                    % Apply first layer matrix
                    L_BC(1:3*sum(Ns(1:(n))), 1:3*sum(Ns(1:(n+1)))) = horzcat(buff_L, buff_R);
                elseif n == numOfLayers % Last layer
                    NL = Ns(n - 1); NR = Ns(n);
                    buff_R = L{1, n};
                    buff_L = zeros(3*NR, 3*NL);
                    % External BC
                    buff_R((2*NR), :) = S{1, n}(NR, :);
                    buff_R((2*NR+1), :) = S{1, n}((2*NR), :);
                    buff_R((3*NR), :) = S{1, n}((3*NR), :);
                    % Displacement BC
                    % (n-1)-th layer
                    buff_L(1, :) = getIBC(Is{1, (n-1)}(NL, :), NL, 1);
                    buff_L(NR, :) = getIBC(Is{1, (n-1)}(NL, :), NL, 2);
                    buff_L((NR+1), :) = getIBC(Is{1, (n-1)}(NL, :), NL, 3);
                    % n-th layer
                    buff_R(1, :) = getIBC(-Is{1, n}(1, :), NR, 1);
                    buff_R(NR, :) = getIBC(-Is{1, n}(1, :), NR, 2);
                    buff_R((NR+1), :) = getIBC(-Is{1, n}(1, :), NR, 3);
                    % Apply last layer matrix
                    L_BC((3*sum(Ns(1:(n - 1)))+1):3*sum(Ns(1:n)), (3*sum(Ns(1:(n - 2)))+1):(3*sum(Ns(1:n)))) = horzcat(buff_L, buff_R);
                else % Middle layers
                    NL = Ns(n-1);
                    NM = Ns(n);
                    NR = Ns(n+1);
                    buff_L = zeros(3*NM, 3*NL);
                    buff_M = L{1, n};
                    buff_R = zeros(3*NM, 3*NR);

                    % Continuity BC
                    % n-th layer
                    buff_M((2*NM), :) = S{1, n}(NM, :);
                    buff_M((2*NM+1), :) = S{1, n}((2*NM), :);
                    buff_M(3*NM, :) = S{1, n}((3*NM), :);
                    % (n+1)-th layer
                    buff_R((2*NM), :) = -S{1, (n+1)}((1), :);
                    buff_R((2*NM+1), :) = -S{1, (n+1)}((NR+1), :);
                    buff_R((3*NM), :) = -S{1, (n+1)}((2*NR+1), :);
                    % Displacement BC
                    % (n-1)-th layer
                    buff_L(1, :) = getIBC(Is{1, (n-1)}((NL), :), NL, 1);
                    buff_L((NM), :) = getIBC(Is{1, (n-1)}((NL), :), NL, 2);
                    buff_L((NM+1), :) = getIBC(Is{1, (n-1)}((NL), :), NL, 3);
                    % n-th layer
                    buff_M(1, :) = getIBC(-Is{1, n}(1, :), NM, 1);
                    buff_M((NM), :) = getIBC(-Is{1, n}(1, :), NM, 2);
                    buff_M((NM+1), :) = getIBC(-Is{1, n}(1, :), NM, 3);
                    % Apply n-th layer matrix
                    L_BC((3*sum(Ns(1:(n - 1)))+1):3*sum(Ns(1:n)), (3*sum(Ns(1:(n-2)))+1):3*sum(Ns(1:(n+1)))) = horzcat(buff_L, buff_M, buff_R);
                end
            end
        case 2 % Modified Quintanilla's scheme of boundary conditions
            for n = 1:numOfLayers
                if n == 1 % First layer
                    NL = Ns(n); NR = Ns(n + 1);
                    buff_L = L{1, n};
                    buff_R = zeros(3*NL, 3*NR);
                    % External BC
                    buff_L(1, :) = S{1, n}(1, :);
                    buff_L((NL+1), :) = S{1, n}((NL+1), :);
                    buff_L((2*NL+1), :) = S{1, n}((2*NL+1), :);
                    % Continuity BC
                    % 1st layer
                    buff_L((NL), :) = S{1, n}(NL, :);
                    buff_L((2*NL), :) = S{1, n}((2*NL), :);
                    buff_L((3*NL), :) = S{1, n}((3*NL), :);
                    % 2nd layer
                    buff_R((NL), :) = -S{1, (n+1)}((1), :);
                    buff_R((2*NL), :) = -S{1, (n+1)}((NR+1), :);
                    buff_R((3*NL), :) = -S{1, (n+1)}((2*NR+1), :);
                    % Apply first layer matrix
                    L_BC(1:3*sum(Ns(1:(n))), 1:3*sum(Ns(1:(n+1)))) = horzcat(buff_L, buff_R);
                elseif n == numOfLayers % Last layer
                    NL = Ns(n - 1); NR = Ns(n);
                    buff_R = L{1, n};
                    buff_L = zeros(3*NR, 3*NL);
                    % External BC
                    buff_R((NR), :) = S{1, n}(NR, :);
                    buff_R((2*NR), :) = S{1, n}((2*NR), :);
                    buff_R((3*NR), :) = S{1, n}((3*NR), :);
                    % Displacement BC
                    % (n-1)-th layer
                    buff_L(1, :) = getIBC(Is{1, (n-1)}(NL, :), NL, 1);
                    buff_L((NR+1), :) = getIBC(Is{1, (n-1)}(NL, :), NL, 2);
                    buff_L((2*NR+1), :) = getIBC(Is{1, (n-1)}(NL, :), NL, 3);
                    % n-th layer
                    buff_R(1, :) = getIBC(-Is{1, n}(1, :), NR, 1);
                    buff_R((NR+1), :) = getIBC(-Is{1, n}(1, :), NR, 2);
                    buff_R((2*NR+1), :) = getIBC(-Is{1, n}(1, :), NR, 3);
                    % Apply last layer matrix
                    L_BC((3*sum(Ns(1:(n - 1)))+1):3*sum(Ns(1:n)), (3*sum(Ns(1:(n - 2)))+1):(3*sum(Ns(1:n)))) = horzcat(buff_L, buff_R);
                else % Middle layers
                    NL = Ns(n-1);
                    NM = Ns(n);
                    NR = Ns(n+1);
                    buff_L = zeros(3*NM, 3*NL);
                    buff_M = L{1, n};
                    buff_R = zeros(3*NM, 3*NR);

                    % Continuity BC
                    % n-th layer
                    buff_M((NM), :) = S{1, n}(NM, :);
                    buff_M((2*NM), :) = S{1, n}((2*NM), :);
                    buff_M(3*NM, :) = S{1, n}((3*NM), :);
                    % (n+1)-th layer
                    buff_R((NM), :) = -S{1, (n+1)}((1), :);
                    buff_R((2*NM), :) = -S{1, (n+1)}((NR+1), :);
                    buff_R((3*NM), :) = -S{1, (n+1)}((2*NR+1), :);
                    % Displacement BC
                    % (n-1)-th layer
                    buff_L(1, :) = getIBC(Is{1, (n-1)}((NL), :), NL, 1);
                    buff_L((NM+1), :) = getIBC(Is{1, (n-1)}((NL), :), NL, 2);
                    buff_L((2*NM+1), :) = getIBC(Is{1, (n-1)}((NL), :), NL, 3);
                    % n-th layer
                    buff_M(1, :) = getIBC(-Is{1, n}(1, :), NM, 1);
                    buff_M((NM+1), :) = getIBC(-Is{1, n}(1, :), NM, 2);
                    buff_M((2*NM+1), :) = getIBC(-Is{1, n}(1, :), NM, 3);
                    % Apply n-th layer matrix
                    L_BC((3*sum(Ns(1:(n - 1)))+1):3*sum(Ns(1:n)), (3*sum(Ns(1:(n-2)))+1):3*sum(Ns(1:(n+1)))) = horzcat(buff_L, buff_M, buff_R);
                end
            end
        case 3 % Original Mekkaoui's scheme (doesn't work)
            for n = 1:numOfLayers
                if n == 1 % First layer
                    NL = Ns(n); NR = Ns(n + 1);
                    buff_L = L{1, n};
                    buff_R = zeros(3*NL, 3*NR);
                    % External BC
                    buff_L(1, :) = S{1, n}(1, :);
                    buff_L((NL), :) = S{1, n}((NL+1), :);
                    buff_L((2*NL), :) = S{1, n}((2*NL+1), :);
                    % Continuity BC
                    % 1st layer
                    buff_L((NL+1), :) = S{1, n}(NL, :);
                    buff_L((2*NL+1), :) = S{1, n}((2*NL), :);
                    % 2nd layer
                    buff_R((NL+1), :) = -S{1, (n+1)}((NR), :);
                    buff_R((2*NL+1), :) = -S{1, (n+1)}((2*NR+1), :);
                    % Displacement BC
                    % 1st layer
                    buff_L((3*NL), :) = getIBC(Is{1, n}(NL, :), NL, 1);
                    % 2тв layer
                    buff_R((3*NL), :) = getIBC(-Is{1, n}(1, :), NR, 1);
                    % Apply first layer matrix
                    L_BC(1:3*sum(Ns(1:(n))), 1:3*sum(Ns(1:(n+1)))) = horzcat(buff_L, buff_R);
                elseif n == numOfLayers % Last layer
                    NL = Ns(n - 1); NR = Ns(n);
                    buff_R = L{1, n};
                    buff_L = zeros(3*NR, 3*NL);
                    % External BC
                    buff_R((NR+1), :) = S{1, n}(NR, :);
                    buff_R((2*NR+1), :) = S{1, n}((3*NR), :);
                    buff_R((3*NR), :) = S{1, n}((3*NR), :);
                    % Continiouty BC
                    % (n-1)-th layer
                    buff_L((1), :) = S{1, (n-1)}((3*NL), :);
                    % n-th layer
                    buff_R((1), :) = -S{1, (n)}((2*NR+1), :);
                    % Displacement BC
                    % (n-1)-th layer
                    buff_L(NR, :) = getIBC(Is{1, (n-1)}(NL, :), NL, 2);
                    buff_L((2*NR), :) = getIBC(Is{1, (n-1)}(NL, :), NL, 3);
                    % n-th layer
                    buff_R(NR, :) = getIBC(-Is{1, n}(1, :), NR, 2);
                    buff_R((2*NR), :) = getIBC(-Is{1, n}(1, :), NR, 3);
                    % Apply last layer matrix
                    L_BC((3*sum(Ns(1:(n - 1)))+1):3*sum(Ns(1:n)), (3*sum(Ns(1:(n - 2)))+1):(3*sum(Ns(1:n)))) = horzcat(buff_L, buff_R);
                else % Middle layers
                    NL = Ns(n-1);
                    NM = Ns(n);
                    NR = Ns(n+1);
                    buff_L = zeros(3*NM, 3*NL);
                    buff_M = L{1, n};
                    buff_R = zeros(3*NM, 3*NR);
                    % Continuity BC
                    % (n-1)-th layer
                    buff_L(1, :) = S{1, (n-1)}((3*NL), :);
                    % n-th layer
                    buff_M(1, :) = -S{1, n}((2*NM+1), :);
                    buff_M((NM+1), :) = S{1, n}((NM), :);
                    buff_M((2*NM+1), :) = S{1, n}((NM), :);
                    % (n+1)-th layer
                    buff_R((NM+1), :) = -S{1, (n+1)}((NR), :);
                    buff_R((2*NM+1), :) = -S{1, (n+1)}((2*NR+1), :);
                    % Displacement BC
                    % (n-1)-th layer
                    buff_L(NM, :) = getIBC(Is{1, (n-1)}((NL), :), NL, 2);
                    buff_L((2*NM), :) = getIBC(Is{1, (n-1)}((NL), :), NL, 3);
                    % n-th layer
                    buff_M(NM, :) = getIBC(-Is{1, n}(1, :), NM, 2);
                    buff_M((2*NM), :) = getIBC(-Is{1, n}(1, :), NM, 3);
                    buff_M((3*NM), :) = getIBC(Is{1, n}((NM), :), NM, 1);
                    % (n+1)-th layer
                    buff_R((3*NR), :) = getIBC(-Is{1, (n+1)}(1, :), NR, 1);
                    % Apply n-th layer matrix
                    L_BC((3*sum(Ns(1:(n - 1)))+1):3*sum(Ns(1:n)), (3*sum(Ns(1:(n-2)))+1):3*sum(Ns(1:(n+1)))) = horzcat(buff_L, buff_M, buff_R);
                end
            end
        case 4 % Modified Mekkaoui's scheme
            for n = 1:numOfLayers
                if n == 1 % First layer
                    NL = Ns(n); NR = Ns(n + 1);
                    buff_L = L{1, n};
                    buff_R = zeros(3*NL, 3*NR);
                    % External BC
                    buff_L(1, :) = S{1, n}(1, :);
                    buff_L((NL), :) = S{1, n}((NL+1), :);
                    buff_L((2*NL), :) = S{1, n}((2*NL+1), :);
                    % Continuity BC
                    % 1st layer
                    buff_L((NL+1), :) = S{1, n}(NL, :);
                    buff_L((2*NL+1), :) = S{1, n}((2*NL), :);
                    % 2nd layer
                    buff_R((NL+1), :) = -S{1, (n+1)}((1), :);
                    buff_R((2*NL+1), :) = -S{1, (n+1)}((NR+1), :);
                    % Displacement BC
                    % 1st layer
                    buff_L((3*NL), :) = getIBC(Is{1, n}(NL, :), NL, 1);
                    % 2тв layer
                    buff_R((3*NL), :) = getIBC(-Is{1, n}(1, :), NR, 1);
                    % Apply first layer matrix
                    L_BC(1:3*sum(Ns(1:(n))), 1:3*sum(Ns(1:(n+1)))) = horzcat(buff_L, buff_R);
                elseif n == numOfLayers % Last layer
                    NL = Ns(n - 1); NR = Ns(n);
                    buff_R = L{1, n};
                    buff_L = zeros(3*NR, 3*NL);
                    % External BC
                    buff_R((NR+1), :) = S{1, n}(NR, :);
                    buff_R((2*NR+1), :) = S{1, n}((2*NR), :);
                    buff_R((3*NR), :) = S{1, n}((3*NR), :);
                    % Continiouty BC
                    % (n-1)-th layer
                    buff_L((1), :) = S{1, (n-1)}((3*NL), :);
                    % n-th layer
                    buff_R((1), :) = -S{1, (n)}((2*NR+1), :);
                    % Displacement BC
                    % (n-1)-th layer
                    buff_L(NR, :) = getIBC(Is{1, (n-1)}(NL, :), NL, 2);
                    buff_L((2*NR), :) = getIBC(Is{1, (n-1)}(NL, :), NL, 3);
                    % n-th layer
                    buff_R(NR, :) = getIBC(-Is{1, n}(1, :), NR, 2);
                    buff_R((2*NR), :) = getIBC(-Is{1, n}(1, :), NR, 3);
                    % Apply last layer matrix
                    L_BC((3*sum(Ns(1:(n - 1)))+1):3*sum(Ns(1:n)), (3*sum(Ns(1:(n - 2)))+1):(3*sum(Ns(1:n)))) = horzcat(buff_L, buff_R);
                else % Middle layers
                    NL = Ns(n-1);
                    NM = Ns(n);
                    NR = Ns(n+1);
                    buff_L = zeros(3*NM, 3*NL);
                    buff_M = L{1, n};
                    buff_R = zeros(3*NM, 3*NR);
                    % Continuity BC
                    % (n-1)-th layer
                    buff_L(1, :) = S{1, (n-1)}((3*NL), :);
                    % n-th layer
                    buff_M(1, :) = -S{1, n}((2*NM+1), :);
                    buff_M((NM+1), :) = S{1, n}((NM), :);
                    buff_M((2*NM+1), :) = S{1, n}((2*NM), :);
                    % (n+1)-th layer
                    buff_R((NM+1), :) = -S{1, (n+1)}((1), :);
                    buff_R((2*NM+1), :) = -S{1, (n+1)}((NR+1), :);
                    % Displacement BC
                    % (n-1)-th layer
                    buff_L(NM, :) = getIBC(Is{1, (n-1)}((NL), :), NL, 2);
                    buff_L((2*NM), :) = getIBC(Is{1, (n-1)}((NL), :), NL, 3);
                    % n-th layer
                    buff_M(NM, :) = getIBC(-Is{1, n}(1, :), NM, 2);
                    buff_M((2*NM), :) = getIBC(-Is{1, n}(1, :), NM, 3);
                    buff_M((3*NM), :) = getIBC(Is{1, n}((NM), :), NM, 1);
                    % (n+1)-th layer
                    buff_R((3*NR), :) = getIBC(-Is{1, (n+1)}(1, :), NR, 1);
                    % Apply n-th layer matrix
                    L_BC((3*sum(Ns(1:(n - 1)))+1):3*sum(Ns(1:n)), (3*sum(Ns(1:(n-2)))+1):3*sum(Ns(1:(n+1)))) = horzcat(buff_L, buff_M, buff_R);
                end
            end
    end

    % Balancing multilayer system
    switch Balance_scheme
        case 0
            L_BC_balanced = L_BC;
            M_BC_balanced = M_BC;
        case 1
            [T, L_BC_balanced] = balance(L_BC);
            M_BC_balanced = T\M_BC*T;
        case 2
            [T1, T2] = balance2(L_BC, M_BC);
            L_BC_balanced = T1*L_BC*T2;
            M_BC_balanced = T1*M_BC*T2;
    end
    
    % Calculation eigenvalues and eigenvectors of equation L*U = w2*M*U
    [U, w2] = eig(L_BC_balanced, M_BC_balanced);
    switch Balance_scheme
        case 0
            % None
        case 1
            % U = T\U*T;
        case 2
            U = T2*U;
    end
    w2 = diag(w2);
    [w2, inds] = sort(w2);
    w = real(sqrt(w2));
    Vp = w./k;
    Vp = Vp';
    fs = w./(2*pi);
    freqs(i, :) = fs/1000;
    velocities(i, :) = Vp;
    
    % Mode separation
    U = real(U);
    U = U(:, inds);
    U1 = []; U2 = []; U3 = [];
    for n = 1:numOfLayers
        U1 = [U1; U((3*sum(Ns(1:(n-1)))+1):(3*sum(Ns(1:n)) - 2*Ns(n)), :)];
        U2 = [U2; U((3*sum(Ns(1:(n-1)))+Ns(n)+1):(3*sum(Ns(1:n)) - Ns(n)), :)];
        U3 = [U3; U((3*sum(Ns(1:(n-1)))+2*Ns(n)+1):(3*sum(Ns(1:n))), :)];
    end
    
    for j = 1:3*sum(Ns)
        % % Check this plot to watch displacements of the interfaces through the thickness.
        % p1 = plot(U1(:, j), 1:sum(Ns));
        % hold on;
        % p2 = plot(U2(:, j), 1:sum(Ns));
        % hold on;
        % p3 = plot(U3(:, j), 1:sum(Ns));
        % % hold on;
        % % plot([-1 1], [Ns(1) Ns(1)], '-k');
        % % hold on;
        % % plot([-1 1], [2*sum(Ns) 2*sum(Ns)], '-k');
        % hold off;
        % % xlim([-1 1]);
        % ylim([1 sum(Ns)])
        % legend([p1, p2, p3], {'U1', 'U2', 'U3'});
        % xlabel('Normalized displacement');
        % ylabel('Points through 3 sets of the collocation points')
        
        % Modes separation
        if Vp(j) ~= 0 && isfinite(Vp(j))
        % Skip values equal to zero or Inf
            % Values of normalized displacements
            U11 = U1(1, j); U1n = U1(end, j);
            U21 = U2(1, j); U2n = U2(end, j);
            U31 = U3(1, j); U3n = U3(end, j);
            if is_orth
                % We have an orthotropic material with propagation angle equal to 0
                % Uncoupled pure Lamb and SH modes
                if (abs(U21) <= z && abs(U2n) <= z) && ~(abs(U21) <= z && abs(U2n) <= z && (abs(U11) <= z) && abs(U1n) <= z && abs(U31) <= z && abs(U3n) <= z) % Pure Lamb wave
                % if (abs(U21) <= z) && (abs(U2n) <= z) % Pure Lamb wave
                    % Case of nulity displacement at the borders of U2 interface
                    if (abs(U11 - U1n) <= z) && (abs(abs(U31) - abs(U3n)) <= z) && (abs(U1n) >= z || (abs(U31 + U3n) <= z))
                        % Symmetric mode:
                        % Interface U1: displacements on the borders are equal;
                        % Interface U3: modulus of displacements on the borbers are equal;
                        % Displacements on the borders of U1 should be greater zero
                        % or displacements on the borders of U3 should be opposite in sign
                        SM(i, j) = Vp(j);
                    else
                        % Otherwise the mode is antisymmetric
                        ASM(i, j) = Vp(j);
                    end
                    Lwave(i, j) = Vp(j);
                elseif ((abs(U11) <= z) && (abs(U1n) <= z) && (abs(U31) <= z) && (abs(U3n) <= z)) % Pure SH mode
                    % Case of nulity displacement at the borders of U1 and U3 interfaces
                    if (abs(U21 - U2n) <= z)
                        % Symmetric mode:
                        % Interface U1: displacements on the borders are equal;
                        SHSM(i, j) = Vp(j);
                    else
                        % Otherwise the mode is antisymmetric
                        SHASM(i, j) = Vp(j);
                    end
                    SHwave(i, j) = Vp(j);
                end                
            end

            % Coupled pure Lamb and SH modes
            z2 = z^2;
            if (U1(1, j)*U1(end, j) >= -z2)  && (U2(1, j)*U2(end, j) >= -z2) && (U3(1, j)*U3(end, j) <= z2)
                % Symmetric mode
                Sym(i, j) = Vp(j);
            else
                % Antisymmetric mode
                Asym(i, j) = Vp(j);
            end
        end
    end
end
timeend = toc(timestart);
disp(['Calculation has taken ', num2str(timeend, '%.3f'), ' seconds']);
delete(f);

%% Plot

if plot_PV
    % Plot all modes without separation
    if plot_PV_all_modes
        figure();
        plot(freqs, velocities, '.r', 'MarkerSize', 2);
        ylim([0 WNLimit]);
        xlim([0 FreqLimit]);
        xlabel('Frequency, kHz');
        ylabel('Phase velocity, m/s');
        title('Lamb and SH modes');
    end
    
    if plot_PV_separated_A_S_modes
        % Plot Symmetric and antisymmetric modes
        figure();
        p1 = plot(freqs, Sym, '.r', 'MarkerSize', 2);
        p1 = p1(1);
        hold on;
        p2 = plot(freqs, Asym, '.b', 'MarkerSize', 2);
        p2 = p2(1);
        hold off;
        ylim([0 WNLimit]);
        xlim([0 FreqLimit]);
        xlabel('Frequency, kHz');
        ylabel('Phase velocity, m/s');
        title('Symmetric and Antisymmetric modes');
        legend([p1, p2], {'S modes', 'AS modes'});
    end

    if is_orth
        if plot_PV_separated_L_SH_modes
            % Plot Lamb and SH modes
            figure();
            p1 = plot(freqs, Lwave, '.r', 'MarkerSize', 2);
            p1 = p1(1);
            hold on;
            p2 = plot(freqs, SHwave, '.b', 'MarkerSize', 2);
            p2 = p2(1);
            hold off;
            ylim([0 WNLimit]);
            xlim([0 FreqLimit]);
            xlabel('Frequency, kHz');
            ylabel('Phase velocity, m/s');
            title('Lamb and SH modes');
            legend([p1, p2], {'Lamb modes', 'SH modes'});
        end
        
        if plot_PV_L_separated_modes
            % Plot Symmetric and Antisymmetric Lamb modes only
            figure();
            p1 = plot(freqs, SM, '.r', 'MarkerSize', 2);
            p1 = p1(1);
            ylim([0 WNLimit]);
            xlim([0 FreqLimit]);
            hold on;
            p2 = plot(freqs, ASM, '.b', 'MarkerSize', 2);
            p2 = p2(1);
            ylim([0 WNLimit]);
            xlim([0 FreqLimit]);
            xlabel('Frequency, kHz');
            ylabel('Phase velocity, m/s');
            title('Lamb symmetric and antisymmetric modes');
            legend([p1, p2], {'Symmetric', 'Antisymmetric'});
        end
        
        if plot_PV_SH_separated_modes
            % Plot Symmetric and Antisymmetric SH modes only
            figure();
            p3 = plot(freqs, SHSM, '.r', 'MarkerSize', 2);
            p3 = p3(1);
            hold on;
            ylim([0 1.2e4]);
            xlim([0 1000]);
            p4 = plot(freqs, SHASM, '.b', 'MarkerSize', 2);
            p4 = p4(1);
            ylim([0 WNLimit]);
            xlim([0 FreqLimit]);
            xlabel('Frequency, kHz');
            ylabel('Phase velocity, m/s');
            title('Shear horizontal symmetric and antisymmetric modes');
            legend([p3, p4], {'Symmetric', 'Antisymmetric'});
        end
        
        if plot_PV_all_modes_separated
            % Plot Symmetric and Antisymmetric, Lamb and SH modes
            figure();
            p1 = plot(freqs, SM, '.r', 'MarkerSize', 2);
            p1 = p1(1);
            hold on;
            p2 = plot(freqs, ASM, '.b', 'MarkerSize', 2);
            p2 = p2(1);
            hold on;
            p3 = plot(freqs, SHSM, '.g', 'MarkerSize', 2);
            p3 = p3(1);
            hold on;
            p4 = plot(freqs, SHASM, '.k', 'MarkerSize', 2);
            p4 = p4(1);
            ylim([0 WNLimit]);
            xlim([0 FreqLimit]);
            xlabel('Frequency, kHz');
            ylabel('Phase velocity, m/s');
            legend([p1, p2, p3, p4], {'S L', 'AS L', 'S SH', 'AS SH'});
        end
    end
end

if plot_GV
    % Plot all modes without separation
    [F, V] = sortVelocities(freqs, velocities);
    GV = calcGV(F, wavenumber);
    figure();
    plot(F, GV, 'b');
    ylim([0 WNLimit]);
    xlim([0 FreqLimit]);
    xlabel('Frequency, kHz');
    ylabel('Group velocity, m/s');
    title('All modes group velocities');

    % % Plot Symmetric and Antisymmetric Lamb modes only
    % [F_S, V_S] = sortVelocities(freqs, SM);
    % GV_S = calcGV(F_S, wavenumber);
    % [F_A, V_A] = sortVelocities(freqs, ASM);
    % GV_A = calcGV(F_A, wavenumber);
    % figure();
    % p1 = plot(F_S, GV_S, 'r');
    % p1 = p1(1);
    % ylim([0 WNLimit]);
    % xlim([0 FreqLimit]);
    % hold on;
    % p2 = plot(F_A, GV_A, '.b');
    % p2 = p2(1);
    % ylim([0 WNLimit]);
    % xlim([0 FreqLimit]);
    % xlabel('Frequency, kHz');
    % ylabel('Phase velocity, m/s');
    % legend([p1, p2], {'SM', 'ASM'});
end

%% Funtions
function c = transformStiffnessMatrix(C, beta)
c = zeros(6,6);
s = sind(beta);
g = cosd(beta);
c(1,1) = C(1,1)*g^4+C(2,2)*s^4+2*(C(1,2)+2*C(6,6))*s^2*g^2; %#ok<*SAGROW>
c(1,2) = (C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2+C(1,2);
c(1,3) = C(1,3)*g^2+C(2,3)*s^2;
c(1,6) = (C(1,2)+2*C(6,6)-C(1,1))*s*g^3+(C(2,2)-C(1,2)-2*C(6,6))*g*s^3;
c(2,2) = C(1,1)*s^4+C(2,2)*g^4+2*(C(1,2)+2*C(6,6))*s^2*g^2;
c(2,3) = C(2,3)*g^2+C(1,3)*s^2;
c(2,6) = (C(1,2)+2*C(6,6)-C(1,1))*g*s^3+(C(2,2)-C(1,2)-2*C(6,6))*s*g^3;
c(3,3) = C(3,3);
c(3,6) = (C(2,3)-C(1,3))*s*g;
c(4,4) = C(4,4)*g^2+C(5,5)*s^2;
c(4,5) = (C(4,4)-C(5,5))*s*g;
c(5,5) = C(5,5)*g^2+C(4,4)*s^2;
c(6,6) = C(6,6)+(C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2;
end

function L = getLmatrix(c, k, I, D, D2)
L12 = c(4,5)*D2 - c(1,6)*I*k^2 + (c(1,4) + c(5,6))*k*1i*D;
L13 = c(3,5)*D2 - c(1,5)*I*k^2 + (c(1,3) + c(5,5))*k*1i*D;
L23 = c(3,4)*D2 - c(5,6)*I*k^2 + (c(3,6) + c(4,5))*k*1i*D;
L11 = c(5,5)*D2 - c(1,1)*I*k^2 + 2*c(1,5)*k*1i*D;
L22 = c(4,4)*D2 - c(6,6)*I*k^2 + 2*c(4,6)*k*1i*D;
L33 = c(3,3)*D2 - c(5,5)*I*k^2 + 2*c(3,5)*k*1i*D;
L21 = L12;
L31 = L13;
L32 = L23;
L = [L11 L12 L13;...
     L21 L22 L23;...
     L31 L32 L33];
end

function S = getSmatrix(c, k, I, D)
S1 = c(1,3)*k*1i*I + c(3,5)*D;
S2 = c(3,6)*k*1i*I + c(3,4)*D;
S3 = c(3,5)*k*1i*I + c(3,3)*D;
S4 = c(1,4)*k*1i*I + c(4,5)*D;
S5 = c(4,6)*k*1i*I + c(4,4)*D;
S6 = c(4,5)*k*1i*I + c(3,4)*D;
S7 = c(1,5)*k*1i*I + c(5,5)*D;
S8 = c(5,6)*k*1i*I + c(4,5)*D;
S9 = c(5,5)*k*1i*I + c(3,5)*D;
S = [S1 S2 S3;...
     S4 S5 S6;...
     S7 S8 S9];
% S = [S4 S5 S6; ...
%      S7 S8 S9; ...
%      S1 S2 S3];
end

function M = getMmatrix(N, p)
M = eye(N*3).*(-p);
end

function I_bc = getIBC(I_n, N, Pos)
I_bc = zeros(1, 3*N);
I_bc(1, (N*(Pos-1) + 1):N*Pos) = I_n;
end

function L = getLBC (L, S, N)
L(1, 1:3*N) = S(1, 1:3*N);
L(N, 1:3*N) = S(N, 1:3*N);
L((N+1), 1:3*N) = S((N+1), 1:3*N);
L(2*N, 1:3*N) = S(2*N, 1:3*N);
L((2*N+1), 1:3*N) = S((2*N+1), 1:3*N);
L(3*N, 1:3*N) = S(3*N, 1:3*N);
end

function M = getMBC(M, N)
zerosN = zeros(1, 3*N);
M(1, :) = zerosN;
M(N, :) = zerosN;
M((N+1), :) = zerosN;
M((2*N), :) = zerosN;
M((2*N+1), :) = zerosN;
M(3*N, :) = zerosN;
end

function is = isOrthotropic(orientation, angle)
is = 1;
if (angle ~= 0 && angle ~= 90)
    is = 0;
end
for i = 1:length(orientation)
    if orientation(i) ~= 90 && orientation(i) ~= 0
        is = 0;
    end
end
end

function isSym = isSymmetric(ps, Cs, Ors, Ths)
isSym = 1;
i_1 = 1;
i_2 = length(ps);
while i_1 < i_2
    if ~(ps(i_1) == ps(i_2) && isequal(Cs(i_1), Cs(i_2)) && Ors(i_1) == Ors(i_2) && Ths(i_1) == Ths(i_2))
        isSym = 0;
    end
    i_1 = i_1 + 1;
    i_2 = i_2 - 1;
end
end

function GV = calcGV(F, k)
w = F.*pi*2*1000;
[NV1, NV2] = size(F);
GV = NaN(size(F));
for i = 1:NV2
    for j = 2:NV1
        if max(isnan(F((j-1):j, i))) == 1
            continue
        end
        GV(j, i) = (w(j, i) - w((j-1), i))/(k(j) - k(j-1));
    end
end
GV(1, :) = GV(2, :);
end

function [F, V] = sortVelocities(F, V)
NV1 = size(V, 1);
NV2 = size(V, 2);
bunch_size = 4;
max_err = 0.2;
for i = 1:NV2
    for j = (bunch_size+1):(NV1-2)
        if max(isnan(V((j-bunch_size):j, i)))
            continue
        elseif max(isnan(V((j-bunch_size):j, i))) == 0
            prev_F = F((j-bunch_size+1):j, i);
            prev_V = V((j-bunch_size+1):j, i);
            next_F = F((j+1):(j+2), :);
            next_V = V((j+1):(j+2), :);
            ind = findProlongation(prev_F, prev_V, next_F, next_V, i, max_err);
            if ind
                buff = F((j+1):end, i);
                F((j+1):end, i) = F((j+1):end, ind);
                F((j+1):end, ind) = buff;
                buff = V((j+1):end, i);
                V((j+1):end, i) = V((j+1):end, ind);
                V((j+1):end, ind) = buff;
                % figure();
                % plot(F, V);
                % legend;
                % ylim([0 8000]);
                % xlim([0 150]);
                % xlabel('Frequency, kHz');
                % ylabel('Phase velocity, m/s');
                % title('All modes');
            end
        end
        
    end
end
end

function ind = findProlongation(prev_F, prev_V, next_F, next_V, init_ind, err)
ind = 0;
max_f_diff = 4*(mean(diff(prev_F)));
f_diff = abs(next_F(1, :) - prev_F(2));
prognosis_1 = NaN(1, size(next_V, 2));

for i = 1:size(next_F, 2)
    if f_diff(i) < max_f_diff && next_F(1, i) > 0 && max(isinf(prev_F)) == 0 && prev_F(end) < next_F(1, i)
        val = interp1(prev_F, prev_V, next_F(1, i), 'spline', 'extrap');
        if val > 0
            prognosis_1(1, i) = val;
        end
    end
end

corrs_1 = abs(prognosis_1 - next_V(1, :))./next_V(1, :);
[~, ind_1] = min(corrs_1);


% err = abs(prognosis_1(ind_1) - next_V(1, ind_1))/next_V(1, ind_1);
if ind_1 ~= init_ind && min(isnan(prognosis_1)) == 0
    prev_V_2 = [prev_V(2:end); next_V(1, ind_1)];
    prev_F_2 = [prev_F(2:end); next_F(1, ind_1)];

    prognosis_2 = interp1(prev_F_2, prev_V_2, next_F(2, ind_1), 'spline', 'extrap');
    corrs_2 = abs(prognosis_2 - next_V(2, ind_1))./next_V(2, ind_1);
    if corrs_2 <= err
        ind = ind_1;
    end
end

end