%% Parameters
N = 64;
N_int = 1;
D = 24;
alpha = 1.3;
sm = 25;
gm = 0.75;
dt = 4e-3;
gamma = 4.258*2*pi;

%% Spiral
[k_sp,T_e] = genvdspiral(N,N_int,D,alpha,sm,gm);
T_e = dt * floor(T_e / dt); % round down T_e to nearest sampling interval
g_sp = @(t) 1/gamma * (k_sp(t) - k_sp(t - dt)) / dt .* (t <= T_e);

%% Gradient ramp down
g0 = 1/gamma * (k_sp(T_e) - k_sp(T_e - dt)) / dt; % determine final gradient value
t_ramp1 = dt * ceil(abs(g0)/sm/dt); % determine minimum ramp time (s.t. slew is not exceeded)
g_rd = @(t) g0 * (1 - t/t_ramp1) .* (t >= 0) .* (t < t_ramp1);

%% Kspace rewinder
k0 = k_sp(T_e - dt) + gamma * sum(g_rd(0:dt:t_ramp1))*dt; % determine total kspace accumulation
[g_kr,t_trap] = trapgrad(abs(k0), sm, gm, dt);
g_kr = @(t) -g_kr(t) * k0/abs(k0);

%% Piecewise gradient formulation
g = @(t) g_sp(t) .* (t >= 0) .* (t <= T_e) + ...
    g_rd(t - T_e) .* (t > T_e) .* (t <= T_e + t_ramp1) + ...
    g_kr(t - T_e - t_ramp1) .* (t > T_e + t_ramp1) .* (t <= T_e + t_ramp1 + t_trap);

t = 0:dt:(T_e + t_ramp1 + 2*t_ramp2 + t_plat + 1);
plot(t,real(g(t)),t,imag(g(t))), hold on
xline(T_e,'--r');
xline(T_e + t_ramp1, '--r');
xline(T_e + t_ramp1 + t_trap, '--r');