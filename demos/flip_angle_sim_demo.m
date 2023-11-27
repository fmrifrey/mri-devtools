%% Add paths
addpath('~djfrey/mri-devtools/matlab/ndsp');
addpath('~djfrey/mri-devtools/matlab/sim');

%% Set parameters
f_mean = 1/100; % mean precession frequency (kHz)
f_std = f_mean / 10; % std of precession frequency (kHz)
N = 1; % number of spins in population
TE = 10; % echo time (ms)
t_end = 20; % total simulation time (ms)
dt = 0.5; % time interval (ms)
dophasecycle = 0; % option to perform phase cycling
varflipfac = 0.33; % variable flip angle factor
opflip = 90; % excitation flip angle (deg)
opflip2 = 120; % refocuser flip angle (deg)
phs_rf1 = 0; % rf1 tx phase (rad)
phs_rf2 = pi/2; % rf2 tx phase (rad)
phs_rx = 0; % rx phase (rad)

%% Set up experiment
t = 0:dt:t_end; % timing array (ms)

% create flip angle schedule
rflocs = [TE/2,TE:TE:t_end];
nechoes = length(rflocs) - 1;
rfangs = zeros(size(rflocs));
rfangs(1) = opflip * exp(1i*phs_rf1);
for echon = 0:nechoes-1
    if (dophasecycle && echon > 1)
        rf2fac = varflipfac + floor(echon/2 - 1) / floor(nechoes/2 - 1) * (1 - varflipfac);
    elseif (~dophasecycle && echon > 0)
        rf2fac = varflipfac + (echon - 1) / (nechoes - 1) * (1 - varflipfac);
    else
        rf2fac = 1;
    end
    if dophasecycle
        rf2fac = rf2fac * (-1)^echon;
    end
    rfangs(echon + 2) = rf2fac * opflip2 * exp(1i*phs_rf2);
end

% create flips array
flips = zeros(size(t));
for i = 1:length(rflocs)
    flips(t == rflocs(i)) = rfangs(i);
end

% Simulate the magnetization
M = flipsim(t,flips,'N',N,'f_mu',f_mean,'f_sd',f_std);

% Plot results
Mpath = [];
for i = 1:length(t)

    subplot(4,2,1:4)
    
    % loop through spins
    for n = 1:N
        % plot the magnetization
        plot3([0,0.75*M(1,i,n)], [0,0.75*M(2,i,n)], [0,0.75*M(3,i,n)], 'r', 'Linewidth', 0.5);
        hold on
    end
    
    % plot the net mag and axis
    M_net = mean(M(:,:,:),3);
    plot3([0,M_net(1,i)], [0,M_net(2,i)], [0,M_net(3,i)], 'r', 'Linewidth', 3);
    
    % plot the path
    if any(abs(flips(i)) > 0, 2)
        for j = 0:20
            R0 = rmat3D('xy',[real(flips(i)), imag(flips(i))]/180*pi);
            Rj = rmat3D('xy',[real(flips(i)), imag(flips(i))]/180*pi*j/20);
            Mpath = [Mpath, Rj/R0 * M_net(:,i)];
        end
    else
        Mpath = [Mpath, M_net(:,i)];
    end
    plot3(Mpath(1,:),Mpath(2,:),Mpath(3,:),'--','Color',[0.5,0.5,0.5]);
    
    % plot the axises
    plot3([-1,1],[0,0],[0,0],'-k');
    text(-1.1,0,0,'-x'); text(1.1,0,0,'x');
    plot3([0,0],[-1,1],[0,0],'-k');
    text(0,-1.1,0,'-y'); text(0,1.1,0,'y');
    plot3([0,0],[0,0],[-1,1],'-k');
    text(0,0,-1.1,'-z'); text(0,0,1.1,'z');
    zlim([-1,1])
    hold off
    axis off
    
    % plot the flip angles
    subplot(4,2,5:6)
    plot(t,real(flips),t,imag(flips)); hold on
    if isvar('xl'), delete(xl), end
    xl = xline(t(i),'--r'); hold off
    title('flip angles')
    legend('x','y','time bar');
    
    % plot the receiever signal
    subplot(4,2,7:8)
    plot(t(1:i), [cos(phs_rx), sin(phs_rx), 0] * M_net(:,1:i)); hold on
    xlim([0,t(end)]); ylim([-1,1]);
    title('Receiever signal')
    
    % display the time
    sgtitle(sprintf('t = %.1fms', t(i)));
    
    % pause and reset
    pause(0.01)
    drawnow
    hold off
    
end