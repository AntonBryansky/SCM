clearvars;
close all;
%% Layer properties
% CFRP laminate properties T300 woven
p = 1437; % (kg/m^3)
S(1,1) = 59.2;	S(1,2) = 8;     S(1,3) = 16.8; % stiffness (GPa)
S(2,1) = 8;	    S(2,2) = 59.2;	S(2,3) = 16.8;
S(3,1) = 2.6;	S(3,2) = 2.6;	S(3,3) = 9.7;
                                                S(4,4) = 8.4;
                                                                S(5,5) = 8.4;
                                                                               S(6,6) = 6.1;
S = S*1e9; % GPa -> Pa
h = 2*1e-3; % Thickness of the layer, m

%% General parameters
% Mechanical properties of layer
density = p; % density, kg/m^3
C = S; % stiffness matrix, Pa
Orientation = 0; % orientation angle of the layer, degree
thickness = h; % thicknesses of the layer, m
propagationAngle = 0; % Wave propagation direction relative to layer orientation angles, degree
% Differential matrices parameters, number of collocation points
N = 11; % Amount of collocation points. Probably, it is good to take odd number to include zero (middle) point in thickness.
% This value is used instead of zero to avoid computational errors
z = 1e-7; % Error of calculation, depends on the number of collocation points

FreqLimit = 1000; % Frequency limit for plots, kHz
% Set range of wavenumbers
WNLimit = 12000;
WNAmount = 5000;
wavenumber = 0:(WNLimit/WNAmount):WNLimit; wavenumber(1) = 1e-6;

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

%% Transform stiffness matrix relatively to the propagation angle
beta = Orientation - propagationAngle; % angle between the main direction and propagation angle, usually equals to zero, degree
c = transformStiffnessMatrix(C, beta); % transformed stiffness matrix

%% Calculation
% Chebyshev differentiation matrices
[x, D] = chebdif(N, 1); % Chebyshev differentiation matrices of 1st order for N amount of collocations points in diapason [-1, ..., 1]
% Calculation new range according to plate thickness
X_min = -thickness/2;
X_max = thickness/2;
X_c = 0.5*(X_min + X_max);
X = 0.5*(X_max - X_min)*x + X_c;
DM = D(:,:,1).*(2/(X_max - X_min)); % Matrices rescaling
DM2 = DM^2; % Chebyshev differentiation matrix of 2nd order
I = eye(N); % Identity matrix

% M matrix is independent of wavenumber
M = getMmatrix(N, density); % It is possible to exclude M matrix calculation from cycle.
M = setMBC(M, N); % Boundary conditions for M matrix

% Output phase velocities
velocities = NaN(length(wavenumber),N*3);
freqs = NaN(length(wavenumber), N*3);

Lwave = NaN(length(wavenumber), N*3);
SHwave = NaN(length(wavenumber), N*3);
Sym = NaN(length(wavenumber), N*3);
Asym = NaN(length(wavenumber), N*3);

SM = NaN(length(wavenumber), N*3);
ASM = NaN(length(wavenumber), N*3);
SHSM = NaN(length(wavenumber), N*3);
SHASM = NaN(length(wavenumber), N*3);

% Calculation dispersion curves
f = waitbar(0, 'Starting calculations...');
timestart = tic;
for i = 1:length(wavenumber)
    progr = i/length(wavenumber);
    waitbar(progr, f, strcat(num2str(round(progr*100, 3)), '%'))
    % calculating L, M and S matrices
    k = wavenumber(i);
    L = getLmatrix(c, k, I, DM, DM2);
    S = getSmatrix(c, k, I, DM);
    
    % Setting boundary and interface conditions
    % Boundary conditions for L matrix
    L = setLBC(L, S, N);
    
    % Calculation eigenvalues and eigenvectors of equation L*U = w2*M*U
    [U, w2] = eig(L, M, 'qz');
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
    U1 = U(1:N, :);
    U2 = U((N+1):2*N, :);
    U3 = U((2*N+1):3*N, :);
    
    for j = 1:3*N
        % Check this plot to see displacements of the interfaces through the thickness.
        % p1 = plot(U1(:, j), 1:N);
        % hold on;
        % p2 = plot(U2(:, j), 1:N);
        % hold on;
        % p3 = plot(U3(:, j), 1:N);
        % % hold on;
        % % plot([-1 1], [N N], '-k');
        % % hold on;
        % % plot([-1 1], [2*N 2*N], '-k');
        % hold off;
        % xlim([-1 1]);
        % ylim([1 N])
        % legend([p1, p2, p3], {'U1', 'U2', 'U3'});
        % xlabel('Normalized displacement');
        % ylabel('Points through 3 sets of the collocation points')
            
        % The value we supposed to be zero or an infitisinal value for approximately equal values
        
        % Modes separation
        if Vp(j) ~= 0 && isfinite(Vp(j)) % Skip values equal to zero or Inf
            % Values of normalized displacements
            U11 = U1(1, j); U1n = U1(end, j);
            U21 = U2(1, j); U2n = U2(end, j);
            U31 = U3(1, j); U3n = U3(end, j);
            if (propagationAngle == 0 || propagationAngle == 90) && (Orientation == 0 || Orientation == 90)
                % Case of uncoupled pure Lamb and SH modes (orthogonal directions)
                if (abs(U21) <= z) && (abs(U2n) <= z) % Pure Lamb wave
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
                elseif (abs(U11) <= z) && (abs(U1n) <= z) && (abs(U31) <= z) && (abs(U3n) <= z) % Pure SH mode
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
            if (U1(1, j)*U1(end, j) >= 0)  && (U2(1, j)*U2(end, j) >= 0) && (U3(1, j)*U3(end, j) < z)
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
        title('Lamb sytmmetric and antisymmetric modes');
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
% if  c(1,6) == 0
%     c(1,6) = 1;
%     c(2,6) = 1;
%     c(3,6) = 1;
%     c(4,5) = 1;
% end
end

function L = getLmatrix(c, k, I, D, D2)
L11 = c(5,5)*D2 - c(1,1)*I*k^2 + 2*c(1,5)*k*1i*D;
L12 = c(4,5)*D2 - c(1,6)*I*k^2 + (c(1,4) + c(5,6))*k*1i*D;
L13 = c(3,5)*D2 - c(1,5)*I*k^2 + (c(1,3) + c(5,5))*k*1i*D;
L21 = L12;
L22 = c(4,4)*D2 - c(6,6)*I*k^2 + 2*c(4,6)*k*1i*D;
L23 = c(3,4)*D2 - c(5,6)*I*k^2 + (c(3,6) + c(4,5))*k*1i*D;
L31 = L13;
L32 = L23;
L33 = c(3,3)*D2 - c(5,5)*I*k^2 + 2*c(3,5)*k*1i*D;
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
% S = [S4 S5 S6;... % S23
%      S7 S8 S9;... % S33
%      S1 S2 S3];   % S13
end

function M = getMmatrix(N, p)
M = eye(N*3).*(-p);
end

function L = setLBC (L, S, N)
L(1, 1:3*N) = S(1, 1:3*N);
L(N, 1:3*N) = S(N, 1:3*N);
L((N+1), 1:3*N) = S((N+1), 1:3*N);
L(2*N, 1:3*N) = S(2*N, 1:3*N);
L((2*N+1), 1:3*N) = S((2*N+1), 1:3*N);
L(3*N, 1:3*N) = S(3*N, 1:3*N);
end

function M = setMBC(M, N)
M(1, 1) = 0;
M(N, N) = 0;
M((N+1), (N+1)) = 0;
M(2*N, 2*N) = 0;
M((2*N+1), (2*N+1)) = 0;
M(3*N, 3*N) = 0;
end

function [x, DM] = chebdif(N, M)
I = eye(N);                          % Identity matrix.
L = logical(I);                      % Logical identity matrix.
n1 = floor(N/2); n2  = ceil(N/2);     % Indices used for flipping trick.
k = [0:N-1]';                        % Compute theta vector.
th = k*pi/(N-1);
x = sin(pi*[N-1:-2:1-N]'/(2*(N-1))); % Compute Chebyshev points.
T = repmat(th/2,1,N);
DX = 2*sin(T'+T).*sin(T'-T);          % Trigonometric identity.
DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];   % Flipping trick.
DX(L) = ones(N,1);                       % Put 1's on the main diagonal of DX.
C = toeplitz((-1).^k);               % C is the matrix with
C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2;     % entries c(k)/c(j)
C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;
Z = 1./DX;                           % Z contains entries 1/(x(k)-x(j))
Z(L) = zeros(N,1);                      % with zeros on the diagonal.
D = eye(N);                          % D contains diff. matrices.

for ell = 1:M
    D = ell*Z.*(C.*repmat(diag(D),1,N) - D); % Off-diagonals
    D(L) = -sum(D');                            % Correct main diagonal of D
    DM(:,:,ell) = D;                                   % Store current D in DM
end
end

function [F, V] = sortVelocities(F, V)
NV1 = size(V, 1);
NV2 = size(V, 2);
bunch_size = 2;
max_err = 0.5;
for i = 1:NV2
    for j = (bunch_size+50):(NV1-2)
        if max(isnan(V((j-bunch_size):j, i)))
            continue
        elseif max(isnan(V((j-bunch_size):j, i))) == 0
            prev_F = F((j-bunch_size+1):j, i);
            prev_V = V((j-bunch_size+1):j, i);
            next_F = F((j+1):(j+2), :);
            next_V = V((j+1):(j+2), :);
            ind = findProlongation(prev_F, prev_V, next_F, next_V, j, max_err);
            if ind
                buff = F((j+1):end, i);
                F((j+1):end, i) = F((j+1):end, ind);
                F((j+1):end, ind) = buff;
                buff = V((j+1):end, i);
                V((j+1):end, i) = V((j+1):end, ind);
                V((j+1):end, ind) = buff;
            end
        end

    end
end
end

function ind = findProlongation(prev_F, prev_V, next_F, next_V, init_ind, err)
ind = 0;
max_f_diff = 3*(diff(prev_F));
f_diff = abs(next_F(1, :) - prev_F(2));
for i = 1:size(next_F, 2)
    if f_diff(i) > max_f_diff
        next_F(:, i) = nan(size(next_F(:, i)));
    end
end

delta_1 = (prev_V(2) - prev_V(1))/(prev_F(2) - prev_F(1));
prognosis_1 = prev_V(2) + (next_F(1, :) - prev_F(2))*delta_1;
corrs_1 = 100*abs(prognosis_1 - next_V(1, :))./next_V(1, :)./next_F(1, :);
[~, ind_1] = min(corrs_1);
if ind_1 ~= init_ind
    delta_2 = (next_V(1, ind_1) - prev_V(2))/(next_F(1, ind_1) - prev_F(2));
    prognosis_2 = next_V(2, ind_1) + (next_F(2, :) - next_F(2, ind_1))*delta_2;
    corrs_2 = 100*abs(prognosis_2 - next_V(2, :))./next_V(2, :);
    [~, ind_2] = min(corrs_2);
    if ind_1 == ind_2 && corrs_2(ind_2) <= err
        ind = ind_2;
    end
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