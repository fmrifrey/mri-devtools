function M = flipsim(t,flip,varargin)

    defaults = struct( ...
        'N', 1, ...
        'f_mu', 1/20, ...
        'f_sd', 1/100, ...
        'T1', 1500, ...
        'T2', 200 ...
        );

    args = vararginparser(defaults, varargin{:});

    if args.N == 1
        args.f_sd = 0;
    end

    % create population of precession frequencies
    fn = args.f_mu + args.f_sd*randn(args.N,1);

    % initialize magnetization
    M = zeros(3, length(t), args.N); % vector component x time x spin index
    M0 = [0;0;1];

    % loop through time points
    for i = 1:length(t)

        % determine the time differential
        if i == 1
            dt = 0;
        else
            dt = t(i) - t(i-1);
        end

        % loop through spins
        for n = 1:args.N

            % determine the previous magnetization
            if i == 1
                Min1 = M0;
            else
                Min1 = M(:,i-1,n);
            end

            % rotation due to free precession
            R_fp = rmat3D('z',2*pi*fn(n)*dt);

            % rotation due to rf flipping
            R_flip = rmat3D('xy',[real(flip(i)),imag(flip(i))] * 1/180*pi);

            % determine magnetization
            M(:,i,n) = R_flip * R_fp * Min1;

        end
    end

end