classdef ellipsoidph
% classdef ellipsoidph()
%
% Part of fmrifrey/mri-devtools software package by David Frey (2023)
%   git@github.com:fmrifrey/mri-devtools.git
%
% Description: Creates an ellipsoid-based digital phantom with analytical
%   k-space signal calculation based on work from Guan Koay, C., Sarlls,
%   J.E., Özarslan, E., (2007) Three-dimensional analytical magnetic
%   resonance imaging phantom in the Fourier domain, Magn. Reson. Med.
%   58(2):430-436, https://doi.org/10.1002/mrm.21292
%
%   ellispoids are stored as 1x10 vectors in the following format:
%   ellipsoidph.e = [A,a,b,c,x0,y0,z0,phi,theta,psi]
%           - 'A':
%               - intensity of ellipsoid
%               - float describing intensity
%           - 'a':
%               - length of x semiaxis
%               - float describing length as fraction of x-fov (range [0,1])
%           - 'b':
%               - length of y semiaxis
%               - float describing length as fraction of y-fov (range [0,1])
%           - 'c':
%               - length of z semiaxis
%               - float describing length as fraction of z-fov (range [0,1])
%           - 'x0':
%               - centroid x coordinate
%               - float describing x coordinate as fraction of x-fov (range
%                   [-1,1])
%           - 'y0':
%               - centroid y coordinate
%               - float describing y coordinate as fraction of y-fov (range
%                   [-1,1])
%           - 'z0':
%               - centroid z coordinate
%               - float describing y coordinate as fraction of z-fov (range
%                   [-1,1])
%           - 'phi':
%               - 1st rotation euler angle (about z)
%               - float describing euler angle in deg
%           - 'theta':
%               - 2nd rotation euler angle (about x)
%               - float describing euler angle in deg
%           - 'psi':
%               - 3rd rotation euler angle (about z)
%               - float describing euler angle in deg
%
%
% Methods:
%
%   function obj = ellipsoidph()
%
%       Description: Constructor function for ellipsoidph class
%
%       Variable input arguments:
%           - 'template':
%               - 3D phantom template name
%               - 'yu-ye-wang', 'shepp-logan', or 'modified-shepp-logan'
%               - no default
%           - 'e':
%               - 'e' matrix for 3D ellipsoids
%               - Nx10 matrix containing ellipsoid parameters in each row
%               - no default
%           - 'fov':
%               - field of view (cm) for kspace sampling calculation
%               - 3-element vector describing fov along x,y,z
%               - default is [16,16,16]
%       Function Output:
%           - obj:
%               - returned output ellipsoidph object
%
%   function [obj,n] = addlayer(obj,im,varargin)
%
%       Description: Function to add a layer to the ellipsoidph class
%
%       Static input arguments:
%           - e:
%               - 'e' matrix for layers to add
%               - see 'e' under constructor function for more info
%               - no default
%
%       Function output:
%           - obj:
%               - returned output ellipsoidph object with edits
%           - n:
%               - layer index that was just added
%               - if not returned, function will print the layer number
%
%   function [obj,n] = editlayer(obj,n,varargin)
%
%       Description: Function to edit a layer in the ellipsoidph class
%
%       Static input arguments:
%           - n:
%               - layer to edit
%               - integer describing layer index
%               - no default, must specify layer
%
%       Variable input arguments:
%           - 'A','a','b','c','x0','y0','z0','phi','theta','psi':
%               - see paramter definitions under ellipsoidph() description
%               - no default (parameter will not be overwritten unless
%                   specified)
%
%       Function output:
%           - obj:
%               - returned output ellipsoidph object with edits
%
%   function obj = swaplayers(obj,n1,n2)
%       
%       Description: Function to swap layers within the ellipsoidph object
%       
%       Static input arguments:
%           - n1:
%               - first layer to swap
%               - integer describing layer index
%               - no default
%           - n2:
%               - second layer to swap
%               - integer describing layer index
%               - no default
%       
%       Function output:
%           - obj:
%               - returned output ellipsoidph object with edits
%
%   function S = getsignal(obj,k)
%   
%       Description: Function to get analytical solution for signal at
%           kspace sample locations
%
%       Static input arguments:
%           - k:
%               - kspace sampling locations
%               - Nx3 vector representing sampling location in [kx,ky,kz]
%               - no defaults
%  
%       Function output:
%           - S:
%               - kspace signal at locations
%               - Nx1 vector representing complex signal
%
%   function im = getimage(obj,dim)
%
%   Description: Function to get discrete 3D image array of phantom
%
%       Static input arguments:
%           - dim:
%               - image dimensions
%               - 3-element vector representing voxel dimensions along
%                   [x,y,z]
%               - no defaults
%  
%       Function output:
%           - im:
%               - 3D image array
%
    
    properties
        e % Ellipsoid properties
        fov % Field of view (cm)
    end
    
    methods

        function obj = ellipsoidph(varargin)
            
            % Set defaults and parse through variable input arguments
            defaults = struct('template', [], ...
                'e', [], ...
                'fov', [16,16,16]);
            args = vararginparser(defaults, varargin{:});

            % Check that e matrix has correct size
            if ~isempty(args.e) && size(args.e,2) ~= 10
                error('e matrix must be Nx10');
            end

            % Create e matrix based on template if specified
            if ~isempty(args.template)
                switch lower(args.template)
                    case 'shepp-logan'
                        args.e = [shepp_logan(); args.e];
                    case 'modified-shepp-logan'
                        args.e = [modified_shepp_logan(); args.e];
                    case 'yu-ye-wang'
                        args.e = [yu_ye_wang(); args.e];
                    otherwise
                        error('invalid template: %s', args.template);
                end
            end

            % Loop through ellipsoids
            for n = 1:size(args.e,1)
                % Add the layer
                obj = obj.addlayer(args.e);
            end

            % Set fov
            obj.fov = args.fov;

        end

        function [obj,n] = addlayer(obj,e)
            obj.e = [obj.e;e]; % append the ellipse matrix
            n = size(obj.e,1);
        end

        function obj = editlayer(obj,n,varargin)

            % Set defaults and parse through variable inputs
            defaults = struct( ...
                'A', obj.e(n,1), ...
                'a', obj.e(n,2), ...
                'b', obj.e(n,3), ...
                'c', obj.e(n,4), ...
                'x0', obj.e(n,5), ...
                'y0', obj.e(n,6), ...
                'z0', obj.e(n,7), ...
                'phi', obj.e(n,8), ...
                'theta', obj.e(n,9), ...
                'psi', obj.e(n,10) ...
                );
            args = vararginparser(defaults, varargin{:});
            
            % Set parameters
            obj.e(n,1) = args.A;
            obj.e(n,2) = args.a;
            obj.e(n,3) = args.b;
            obj.e(n,4) = args.c;
            obj.e(n,5) = args.x0;
            obj.e(n,6) = args.y0;
            obj.e(n,7) = args.z0;
            obj.e(n,8) = args.phi;
            obj.e(n,9) = args.theta;
            obj.e(n,10) = args.psi;

        end

        function obj = swaplayers(obj,n1,n2)

            % Swap the indicies of the E array
            tmp = obj.e(n1,:);
            obj.e(n1,:) = obj.e(n2,:);
            obj.e(n2,:) = tmp;

        end

        function S = getsignal(obj,k)
            
            % Define bessel functions
            Jn = @(n,z) sqrt(pi./(2*z)).*besselj(n+0.5,z); % Spherical bessel function of the nth kind
            Gn = @(n,K) Jn(n,abs(K))./abs(K);
            
            % Initialize the signal vector
            S = zeros(1, size(k,2));
            
            % Loop through ellipsoids
            for n = 1:size(obj.e,1)

                % Get ellipsoid properties
                rho = obj.e(n,1);
                D = diag(obj.fov(:)'/2 .* obj.e(n,2:4));
                rc = obj.fov(:)/2 .* obj.e(n,5:7)';
                R = rmat3D('zxz', pi/180*obj.e(n,10:-1:8));

                % Add fourier transform of the current ellipsoid
                S = S + ...
                    rho*2*pi*abs(det(D))*exp(-1i*2*pi*rc'*k).*Gn(1,2*pi*vecnorm(D*R'*k,2,1));

            end
        end

        function im = getimage(obj,dim)

            % Initialize image array
            im = zeros(dim);

            % Make image grid
            [X,Y,Z] = imgrid(obj.fov, dim);
            r = [X(:),Y(:),Z(:)]';

            % Loop through ellipsoids
            for n = 1:size(obj.e,1)

                % Get properties
                rho = obj.e(n,1);
                D = diag(obj.fov(:)'/2 .* obj.e(n,2:4));
                rc = obj.fov(:)/2 .* obj.e(n,5:7)';
                R = rmat3D('zxz', pi/180*obj.e(n,10:-1:8));

                % Determine region and add amplitude
                ROI = vecnorm(D^-1*R'*(r - rc), 2, 1) <= 1;
                im(ROI) = im(ROI) + rho;

            end

        end

    end
end


%% Phantom templates
function e = shepp_logan

e = modified_shepp_logan;
e(:,1) = [1 -.98 -.02 -.02 .01 .01 .01 .01 .01 .01];

end
      
function e = modified_shepp_logan
%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%      
%         A      a     b     c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e =    [  1  .6900  .920  .810      0       0       0      0      0      0
        -.8  .6624  .874  .780      0  -.0184       0      0      0      0
        -.2  .1100  .310  .220    .22       0       0    -18      0     10
        -.2  .1600  .410  .280   -.22       0       0     18      0     10
         .1  .2100  .250  .410      0     .35    -.15      0      0      0
         .1  .0460  .046  .050      0      .1     .25      0      0      0
         .1  .0460  .046  .050      0     -.1     .25      0      0      0
         .1  .0460  .023  .050   -.08   -.605       0      0      0      0
         .1  .0230  .023  .020      0   -.606       0      0      0      0
         .1  .0230  .046  .020    .06   -.605       0      0      0      0 ];
       
end       

function e = yu_ye_wang
%
%   Yu H, Ye Y, Wang G, Katsevich-Type Algorithms for Variable Radius Spiral Cone-Beam CT
%      
%         A      a     b     c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e =    [  1  .6900  .920  .900      0       0       0      0      0      0
        -.8  .6624  .874  .880      0       0       0      0      0      0
        -.2  .4100  .160  .210   -.22       0    -.25    108      0      0
        -.2  .3100  .110  .220    .22       0    -.25     72      0      0
         .2  .2100  .250  .500      0     .35    -.25      0      0      0
         .2  .0460  .046  .046      0      .1    -.25      0      0      0
         .1  .0460  .023  .020   -.08    -.65    -.25      0      0      0
         .1  .0460  .023  .020    .06    -.65    -.25     90      0      0
         .2  .0560  .040  .100    .06   -.105    .625     90      0      0
        -.2  .0560  .056  .100      0    .100    .625      0      0      0 ];
       
end