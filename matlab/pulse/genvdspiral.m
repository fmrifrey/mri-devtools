function [k,T_e] = genvdspiral(N, N_int, D, alpha, sm, gm)
% Kim DH, Adalsteinsson E, Spielman DM. Simple analytic variable density
% spiral design. Magn Reson Med. 2003;50(1):214–219.

%% Calculate constants
gamma = 4.258*2*pi; % gyromagnetic ratio (rad/G/ms)
n = (1 - (1 - 2/N) ^ (1/alpha)) ^ (-1) / N_int; % # of turns in kspace
omega = 2 * pi * n; % total angular displacement
lambda = N / (2*D); % kspace maximum

%% Amplitude limited regime
% Eq. [4]: calculation of tau
tau_a = @(t) ((gamma * gm) / (lambda * omega) * (alpha + 1) * t) .^ (1/(alpha + 1));

% Eq. [5]: calculation of trajectory end time
T_ea = ((gamma * gm) / (lambda * omega) * (alpha + 1)) ^ (-1);

%% Slew limited regime
% Eq. [7]: calculation of tau
tau_s = @(t) (sqrt((sm * gamma) / (lambda * omega^2)) * (alpha/2 + 1) * t) .^ (1/(alpha/2 + 1));

% Eq. [8]: calculation of trajectory end time
T_es = ((alpha/2 + 1) * sqrt((sm * gamma) / (lambda * omega^2))) ^ (-1);

%% Slew rate and amplitude limited regime
% Eq. [9]: calculation of transition time
T_s2a = ((gm * gamma) / ...
    (lambda * omega * 1/(alpha + 1) * ...
    ((alpha/2 + 1) * sqrt((sm * gamma) / (lambda * omega^2))) ...
    ^ ((alpha + 1)/(alpha/2 + 1)))) ...
    ^ ((alpha + 2)/alpha);

% Determine end time
if T_s2a > T_es
    T_e = T_es;
else
    T_e = T_ea;
end

% Eq. [11]: calculation of tau
tau = @(t) tau_s(t) .* (t >= 0) .* (t <= min(T_s2a, T_es)) + ...
    tau_a(t) .* (t >= T_s2a) .* (t <= T_ea);

% Eq. [2]: calculation of k-space trajectory before ramps
k = @(t) lambda * tau(t).^alpha .* exp(1i*omega * tau(t));

end
