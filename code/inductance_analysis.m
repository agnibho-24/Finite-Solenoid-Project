% inductance_analysis.m
% This script computes and plots the total inductance of a finite solenoid.
% It studies how inductance changes with:
% 1. number of turns N
% 2. pitch
% 3. loop radius a
% 4. wire radius rw
% The script uses self-inductance of each loop plus numerical mutual
% inductance between loop pairs.
clear; clc; close all;

%% Constants
mu0 = 4*pi*1e-7;

%% Base parameters (used unless varied)
a_base   = 0.05;     % loop radius [m]
rw_base  = 1e-3;     % wire radius [m]
ell_base = 0.20;     % total length [m]
N_base   = 60;       % number of turns

%% ================= 1) L vs N =================
Nvals = 10:10:120;
L_N = zeros(size(Nvals));

for k = 1:length(Nvals)
    N = Nvals(k);
    L_N(k) = computeSolenoidInductance(a_base, ell_base, N, rw_base, mu0);
end

figure;
plot(Nvals, 1e3*L_N, 'o-', 'LineWidth', 1.5);
xlabel('Number of turns N');
ylabel('Inductance, L_T (mH)');
title('L_T vs N');
grid on;

%% ================= 2) L vs Pitch =================
pitchVals = linspace(0.002, 0.2, 50);  % avoid overlap (< 2*rw)
L_pitch = zeros(size(pitchVals));
L_asymptote = N_base * mu0 * a_base * (log(8*a_base/rw_base) - 1.75);

for k = 1:length(pitchVals)
    pitch = pitchVals(k);
    L_pitch(k) = computeSolenoidInductancePitch(a_base, pitch, N_base, rw_base, mu0);
end

figure;
plot(1e3*pitchVals, 1e3*L_pitch, 'o-', 'LineWidth', 1.5);
hold on;

% Asymptote line
yline(1e3*L_asymptote, '--r', 'LineWidth', 1.5);

xlabel('Pitch (mm)');
ylabel('Inductance, L_T (mH)');
title(['L_t vs Pitch (N = ', num2str(N_base), ')']);
legend('Calculated L_T', 'Asymptote: N L_c', 'Location', 'northeast');
grid on;

%% ================= 3) L vs Radius a =================
aVals = linspace(0.02, 0.10, 25);
L_a = zeros(size(aVals));

for k = 1:length(aVals)
    a = aVals(k);
    L_a(k) = computeSolenoidInductance(a, ell_base, N_base, rw_base, mu0);
end

figure;
plot(1e3*aVals, 1e3*L_a, 'o-', 'LineWidth', 1.5);
xlabel('Loop radius a (mm)');
ylabel('Inductance, L_T (mH)');
title(['L_T vs Radius a (N = ', num2str(N_base), ')']);
grid on;

%% ================= 4) L vs Wire Radius rw =================
rwVals = linspace(0.0005, 0.005, 25);
L_rw = zeros(size(rwVals));

for k = 1:length(rwVals)
    rw = rwVals(k);
    L_rw(k) = computeSolenoidInductance(a_base, ell_base, N_base, rw, mu0);
end

figure;
plot(1e3*rwVals, 1e3*L_rw, 'o-', 'LineWidth', 1.5);
xlabel('Wire radius (thickness)  r_w (mm)');
ylabel('Inductance, L_T (mH)');
title(['L_T vs Wire Radius (N = ', num2str(N_base), ')']);
grid on;

%% ================= FUNCTIONS =================

function Ltotal = computeSolenoidInductance(a, ell, N, rw, mu0)

    pitch = ell / N;
    z = linspace(-ell/2 + pitch/2, ell/2 - pitch/2, N);

    % Self term
    Lself_one = mu0 * a * (log(8*a/rw) - 1.75);
    Lself_total = N * Lself_one;

    % Mutual term
    Msum = 0;
    for i = 1:N
        for j = i+1:N
            dz = abs(z(i) - z(j));
            Msum = Msum + mutualLoopNumerical(a, dz, mu0);
        end
    end

    Ltotal = Lself_total + 2*Msum;
end

function Ltotal = computeSolenoidInductancePitch(a, pitch, N, rw, mu0)

    ell = N * pitch;
    z = linspace(-ell/2 + pitch/2, ell/2 - pitch/2, N);

    % Self term
    Lself_one = mu0 * a * (log(8*a/rw) - 1.75);
    Lself_total = N * Lself_one;

    % Mutual term
    Msum = 0;
    for i = 1:N
        for j = i+1:N
            dz = abs(z(i) - z(j));
            Msum = Msum + mutualLoopNumerical(a, dz, mu0);
        end
    end

    Ltotal = Lself_total + 2*Msum;
end

function M = mutualLoopNumerical(a, dz, mu0)

    integrand = @(phi) cos(phi) ./ sqrt(2*a^2*(1 - cos(phi)) + dz.^2);

    M = (mu0 * a^2 / 2) * integral(integrand, 0, 2*pi, ...
        'RelTol', 1e-8, 'AbsTol', 1e-11);
end
