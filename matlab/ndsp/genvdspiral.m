function k = genvdspiral(N,N_int,D,r,alpha,smax,gmax)
% Chang C, Glover GH. Variable-density spiral-in/out functional magnetic
% resonance imaging. Magn Reson Med. 2011 May; 65(5):1287-96.
% doi: 10.1002/mrm.22722. Epub 2010 Nov 30. PMID: 21500257;
% PMCID: PMC3126879.

% smax is in units of G/cm/s
% gmax is in units of G/cm

dt = 4e-6; % s
t = 0:dt:1; % s

% Calculate trajectory constants
r0 = r - (1/alpha)^(alpha/(alpha-1));
tau0 = r - (1/alpha)^(1/(alpha-1));
taumax = tau0 + (1 - r0)^(1/alpha);

% Calculate the archimedian spiral
tauA = genspiralA(N_int, N, D, smax, gmax);

% Calculate the variable density spiral
tauV = genspiralV(N_int, N, D, smax, gmax); 

end

function tau = genspiralA(N_int, N, D, smax, gmax)
% Glover GH. Simple analytic spiral K-space algorithm. Magn Reson Med.
% 1999;42(2):412–415.

% Define/calculate constants
gamma = 26754; % rad/s/G
lambda = N_int/D; % 1/cm
beta = smax * gamma / lambda;
a2 = (9/4 * beta)^(1/3);
q = 1.5; % slew overshoot factor

% Calculate approximate slew-limited to amp-limited transition time from
% Eq. [10]
t_s2a = (gmax * 3*gamma/(2*lambda) / a2^2)^3;

% Calculate slew rate limited formulation of tau from Eq. [7]
tau_s = @(t) 1/2 * beta*t.^2 ./ (1/q + beta/(2*a2) .* t.^(4/3)) / (pi*N);

% Calculate gradient amplitude limited formulation of tau from Eq. [12] 
tau_g = @(t) sqrt((tau_s(t_s2a)*pi*N)^2 + 2*gamma/lambda*gmax*(t - t_s2a)) / (pi*N);

% Calculate total readout time for slew-limited regime from Eq. [8]
t_end_s = 2/(3*N_int) * sqrt((N*pi)^3 / (gamma*D*smax));

% Calculate total readout time for amp-limited regime from Eq. [13]
t_end_g = t_s2a + lambda/(2*gamma*gmax) * ((pi*N/N_int)^2 - (tau_s(t_s2a)*pi*N)^2);

% Define final piecewise tau function
tau = @(t) tau_s(t) .* (t <= min(t_s2a, t_end_s)) + ...
    tau_g(t) .* (t > t_s2a) .* (t <= t_end_g);

end

function tau = genspiralV(N_int, N, D, smax, gmax)
% Kim DH, Adalsteinsson E, Spielman DM. Simple analytic variable density
% spiral design. Magn Reson Med. 2003;50(1):214–219.

% Define/calculate constants
gamma = 26754; % rad/s/G
lambda = N_int/D; % 1/cm
alpha = 2; % undersampling factor
omega = 2*pi*N;
a1 = alpha + 1;
a2 = alpha/2 + 1;
c1 = gmax*gamma / (lambda*omega);
c2 = sqrt(smax*gamma / (lambda*omega^2));

% Calculate approximate slew-limited to amp-limited transition time from
% Eq. [9]
t_s2a = (c1 * a1 / (a2 * c2)^(a1/a2)) ^ ((alpha + 2)/alpha);

% Calculate slew rate limited formulation of tau from Eq. [7]
tau_s = @(t) (c2 * a2 * t) .^ (1/a2);

% Calculate gradient amplitude limited formulation of tau from Eq. [4] 
tau_g = @(t) (c1 * a1 * t) .^ (1/a1);

% Calculate total readout time for slew-limited regime from Eq. [8]
t_end_s = 1/(a2 * c2);

% Calculate total readout time for amp-limited regime from Eq. [13]
t_end_g = 1/(a1 * c1);

% Define final piecewise tau function
tau = @(t) tau_s(t) .* (t <= min(t_s2a, t_end_s)) + ...
    tau_g(t) .* (t > t_s2a) .* (t <= t_end_g);

end
