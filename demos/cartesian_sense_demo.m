%% Make a phantom
D = 24; % 3D FOV (cm)
N = 32; % 3D matrix size

e = ellipsoidph('fov',D*ones(1,3),'template','modified-shepp-logan');
im_truth = e.getimage(N*ones(1,3));

%% Simulate fully-sampled data
% Create fully sampled cartesian kspace grid
[Kx_fs,Ky_fs,Kz_fs] = imgrid(N/D*ones(1,3),N*ones(1,3));

% Simulate the signal
S = e.getsignal([Kx_fs(:),Ky_fs(:),Kz_fs(:)]');
kdat = reshape(S(1,:),N*ones(1,3));

% Inverse fft to reconstruct the image
im_rec = ifftc(kdat,1:3);

cfigopen('Full cartesian sampling');

subplot(2,1,1)
orthoview(im_truth)
title('ground truth')

subplot(2,1,2)
orthoview(abs(im_rec))
title('reconstructed image')

%% Simulate undersampling along z (NOT YET WORKING)
usfac = 4; % N/usfac must be an integer
nc = 4; % number of coils
e.c = simsmaps(D,nc,D/2);

% Create cartesian kspace grid
[Kx_us,Ky_us,Kz_us] = imgrid(N/D*ones(1,3),N*[1,1,1/usfac]);

% Simulate the signal
t = tic();
S = e.getsignal([Kx_us(:),Ky_us(:),Kz_us(:)]');
fprintf('time elapsed = %.3fs\n', toc(t));
kdat = reshape(S,[nc,N*[1,1,1/usfac]]);
kdat = permute(kdat,[2:4,1]);

% Inverse fft to reconstruct the image
im_rec = ifftc(kdat,1:2);

cfigopen(sprintf('%dx Cartesian undersampling',usfac))

subplot(nc+1,1,1)
orthoview(im_truth)
title('ground truth')

for i = 1:nc
    subplot(nc+1,1,1+i)
    orthoview(im_rec,'frame',i)
    title(sprintf('coil %d reconstructed image',i))
end

function c = simsmaps(D,nc,w)

PHI = (3 - sqrt(5)) / 2; % 1D golden ratio
    
% Set default sensitivity width
if nargin < 3 || isempty(w)
    w = D/4;
end

c = [ones(nc,1), 3/4*D*[cos(2*pi*(1:nc)'/PHI)/2, ...
    sin(2*pi*(1:nc)'/PHI)/2, ...
    (1:nc)'/nc-1/2], w*ones(nc,1)];

end