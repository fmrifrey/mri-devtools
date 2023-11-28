%% Create the reference images using phantoms
dim = 64*ones(1,3);

% Structural image
tmp = ellipsoidph;
tmp = tmp.addlayer('rho', 1.5, 'd', [0,0,0], 'r', [0.5,0.75,0.5], ...
    'theta', [0,0,0]);
tmp = tmp.addlayer('rho', -0.5, 'd', [0,0,0], 'r', [0.45,0.7,0.45], ...
    'theta', [0,0,0]);
anotomical = tmp.getimage(dim);

% Condition 1 (a little bit of visual and motor activation)
tmp = ellipsoidph;
tmp = tmp.addlayer('rho', 4, 'd', [0,-0.5,-0.25], 'r', [0.2,0.1,0.1], ...
    'theta', [0,0,0]);
tmp = tmp.addlayer('rho', 2.5, 'd', [0.25,0,0.5], 'r', [0.1,0.2,0.1], ...
    'theta', [0,0,0]);
condition1 = tmp.getimage(dim);

% Condition 2 (a little bit of frontal lobe)
tmp = ellipsoidph;
tmp = tmp.addlayer('rho', 3, 'd', [-0.1,0.35,0], 'r', [0.12,0.14,0.1], ...
    'theta', [0,0,25]);
condition2 = tmp.getimage(dim);

%% Make the figure:
% Create the overlayview object
o = overlayview;

% Add layers with specific properties:
o = o.addlayer(anotomical, 'caxis', [0 1.5]);
o = o.addlayer(condition1, 'caxis', [2 5], 'cmap', hot(128));
o = o.addlayer(condition2, 'caxis', [2 5], 'cmap', cool(128));

% Swap the 2nd and 3rd layer (condition 1 and 2)
o = o.swaplayers(2,3);

% Add labels to the layers
o = o.editlayer(2, 'name', 'condition 2');
o = o.editlayer(3, 'name', 'condition 1');

% Show the overlaid images
o.show();