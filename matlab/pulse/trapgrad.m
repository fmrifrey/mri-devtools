function [g, T_e] = trapgrad(dk, sm, gm, dt)

    gamma = 4.258*2*pi; % gyromagnetic ratio (rad/G/ms)

    % calculate area under waveform
    A = abs(dk)/gamma;
    
    % calculate height of waveform
    h = sqrt(A * sm) * dk/abs(dk);
    if h > gm % limit the amplitude to the max gradient
        h = gm;
    end
    
    % calculate the ramp and plateau time and round to nearest sampling int
    t_ramp = dt * ceil(h/sm/dt);
    t_plateau = dt * ceil((A/h - t_ramp)/dt);
    
    % correct the height for rounding errors
    h = A / (t_ramp + t_plateau);
    
    % formulate the piecewise gradient waveform
    t1 = t_ramp;
    t2 = t_ramp + t_plateau;
    t3 = 2*t_ramp + t_plateau;
    g = @(t) h .* (t/t_ramp) .* (t >= 0) .* (t < t1) ...
        + h .* (t >= t1) .* (t < t2) ...
        + h .* (1 - (t - t2)/t_ramp) .* (t >= t2) .* (t <= t3);
    T_e = t3;
    
end