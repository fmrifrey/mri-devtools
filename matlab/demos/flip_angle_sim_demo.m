f_mean = 1/20; % mean precession frequency (kHz)
N = 2; % number of spins in population
if (N == 1)
    f_std = 0; % std of precession frequency (kHz)
else
    f_std = f_mean / 10; % std of precession frequency (kHz)
end
TE = 10; % echo time (ms)
t_end = 200; % ms
t = 0:0.5:t_end; % timing (ms)

% create flip angle schedule
flips = zeros(length(t),2);
flips(t == TE/2,1) = 90; % add a tipdown at time 5
tmp = 1;
for t_tmp = TE:TE:(t_end - mod(t_end,TE)) % at times 10, 20, 30, 40
    flips(t == t_tmp,2) = tmp*120; % add a refocuser
    tmp = tmp*-1;
end

%% Run sim
% create population of precession frequencies
fn = f_mean + f_std*randn(N,1);

% initialize magnetization
M = zeros(3,length(t),N); % vector component x time x spin index
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
    for n = 1:N
        
        % determine the previous magnetization
        if i == 1
            Min1 = M0;
        else
            Min1 = M(:,i-1,n);
        end
        
        % rotation due to free precession
        R_fp = rmat3D('z',2*pi*fn(n)*dt);
        
        % rotation due to rf flipping
        R_flip = rmat3D('xy',flips(i,:)/180*pi);
        
        % determine magnetization
        M(:,i,n) = R_flip * R_fp * Min1;
        
    end
end

%% Plot results
path = [];
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
    if any(abs(flips(i,:)) > 0, 2)
        for j = 0:20
            path = [path, rmat3D('xy',flips(i,:)/180*pi*j/20)*inv(rmat3D('xy',flips(i,:)/180*pi))*M_net(:,i)];
        end
    else
        path = [path, M_net(:,i)];
    end
    plot3(path(1,:),path(2,:),path(3,:),'--','Color',[0.5,0.5,0.5]);
    
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
    plot(t,flips); hold on
    if isvar('xl'), delete(xl), end
    xl = xline(t(i),'--r'); hold off
    yticks([90 180])
    title('flip angles')
    legend('x','y','time');
    
    % plot the net magnetization
    subplot(4,2,7:8)
    plot(t(1:i),M_net(1,1:i,:)); hold on
    plot(t(1:i),M_net(2,1:i,:)); hold off
    xlim([0,t(end)]); ylim([-1,1]);
    title('net magnetization')
    legend('x','y');
    
    % display the time
    sgtitle(sprintf('t = %.1fms', t(i)));
    
    % pause and reset
    pause(0.1)
    hold off
    
end