% Set fov and dim
fov = ones(1,3);
dim = 64*ones(1,3);

% Create phantom from template
ph = ellipsoidph('fov',fov,'template','modified-shepp-logan');

% Add an extra ellipsoid to the phantom
ph = ph.addlayer('rho',0.05,'d',[0,0,0],'r',[0.1,0.1,0.4],'theta',[0,0,45]);

% Get the image
im = ph.getimage(dim);

% Set up kspace sampling
[Kx,Ky,Kz] = imgrid(dim./fov,dim);
k = [Kx(:),Ky(:),Kz(:)];
clear Kx Ky Kz

% Get the kspace signal
S = ph.getsignal(k);

% Grid and reconstruct the data
S_cart = reshape(S,dim);
im_recon = ifftc(S_cart);

% Plot results
subplot(2,1,1)
lbview(im);
title('reference');
subplot(2,1,2)
lbview(im_recon)
title('reconstructed');