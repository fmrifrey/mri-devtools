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
%           - im:
%               - image to display
%               - a float/double 3D image array or name of a .nii file
%               - only first frame will be used if there is more than 1
%               - no default; required argument
%
%       Variable input arguments:
%           - 'e':
%               - 'e' matrix for layers to add
%               - see 'e' under constructor function for more info
%               - no default
%           - 'rho':
%               - magnitude/proton density of ellipsoid layer
%               - scalar value describing the intensity
%               - no default
%           - 'r':
%               - radii of ellipsoid layer
%               - 3-element vector describing the ellipsoid radii as a
%                   fraction of fov
%               - no default
%           - 'd':
%               - central displacement of ellipsoid layer
%               - 3-element vector describing the ellipsoid offset as a
%                   fraction of fov
%               - no default
%           - 'theta':
%               - rotations of ellipsoid layer
%               - 3-element vector describing the rotations about the
%                   origin in degrees
%               - no default
%
%       Function output:
%           - obj:
%               - returned output imoverlay object with edits
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
%               - no default
%
%       Variable input arguments:
%           - Can change all layer properties (same arguments as in
%               addlayer(), other than 'e')
%
%       Function output:
%           - obj:
%               - returned output imoverlay object with edits
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
%   function e = getmatrix(obj)
%       
%       Description: Function that returns the ellipsoid matrix
%           representation of the phantom
%
%       Function output:
%           - e:
%               - ellipsoid matrix
%               - see 'e' under constructor for more info
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
        E % Ellipsoid properties
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
                [obj,~] = obj.addlayer('e', args.e);
            end

            % Set fov
            obj.fov = args.fov;

        end

        function [obj,n] = addlayer(obj,varargin)

            % Set defaults and parse through variable input arguments
            defaults = struct('e', [], ...
                'rho', [], ...
                'r', [], ...
                'd', [], ...
                'theta', []);
            args = vararginparser(defaults, varargin{:});

            % Get layer index
            n = length(obj.E) + 1;

            % Initialize ellipsoid structure
            obj.E{n} = struct( ...
                'rho', [], ... % Amplitude (proton density)
                'r', [], ... % Radii (x,y,z) (fraction of fov/2)
                'd', [], ... % Central displacement (x,y,z) (fraction of fov/2)
                'theta', []); % Rotation angles (x,y,z) (radians)

            % Update parameters
            if ~isempty(args.rho)
                obj = obj.editlayer(n,'rho',args.rho);
            elseif ~isempty(args.e)
                obj = obj.editlayer(n,'rho',args.e(n,1));
            end

            if ~isempty(args.r)
                obj = obj.editlayer(n,'r',args.r);
            elseif ~isempty(args.e)
                obj = obj.editlayer(n,'r',args.e(n,2:4));
            end

            if ~isempty(args.d)
                obj = obj.editlayer(n,'d',args.d);
            elseif ~isempty(args.e)
                obj = obj.editlayer(n,'d',args.e(n,5:7));
            end

            if ~isempty(args.theta)
                obj = obj.editlayer(n,'theta',args.theta);
            elseif ~isempty(args.e)
                obj = obj.editlayer(n,'theta',args.e(n,8:10));
            end

        end

        function obj = editlayer(obj,n,varargin)

            % Set defaults and parse through variable inputs
            defaults = obj.E{n};
            args = vararginparser(defaults, varargin{:});
            
            % Set parameters
            obj.E{n}.rho = args.rho;
            obj.E{n}.r = args.r;
            obj.E{n}.d = args.d;
            obj.E{n}.theta = args.theta;

        end

        function obj = swaplayers(obj,n1,n2)

            % Swap the indicies of the E array
            tmp = obj.E{n1};
            obj.E{n1} = obj.E{n2};
            obj.E{n2} = tmp;

        end

        function e = getmatrix(obj)
            
            % Initialize e matrix
            e = zeros(length(obj.E),10);

            % Loop through ellipsoids
            for n = 1:length(obj.E)
                % Set the row values
                e(n,:) = [obj.E{n}.rho, ...
                    obj.E{n}.r, ...
                    obj.E{n}.d, ...
                    obj.E{n}.theta];
            end

        end

        function S = getsignal(obj,k)

            % Initialize signal vector
            S = zeros(size(k,1),1);

            % Loop through ellipsoids
            for n = 1:length(obj.E)
                
                % Get indexed ellipsoid
                En = obj.E{n};

                % Get properties
                rho = En.rho;
                abc = obj.fov(:)'/2 .* En.r;
                d = obj.fov(:)'/2 .* En.d;
                A = rmat3D('zyx', pi/180*En.theta(end:-1:1));

                % Calculate kspace signal
                kn0 = k(vecnorm(k,2,2)~=0,:);
                ktd = kn0*d';
                K = vecnorm(kn0*A'.*abc,2,2);
                S(vecnorm(k,2,2)~=0) = S(vecnorm(k,2,2)~=0) + ...
                    rho * abs(det(A')) * prod(abc) * ...
                    exp(-1i * 2*pi*ktd) .* ...
                    ( sin(2*pi*K) - 2*pi*K.*cos(2*pi*K) ) ./ (2*pi^2*K.^3);
                S(vecnorm(k,2,2)==0) = S(vecnorm(k,2,2)==0) + ...
                    rho*4/3*pi*prod(abc);

            end

        end

        function im = getimage(obj,dim)

            % Initialize image array
            im = zeros(dim);

            % Make image grid
            [X,Y,Z] = imgrid(obj.fov, dim);
            x = [X(:),Y(:),Z(:)];

            % Loop through ellipsoids
            for n = 1:length(obj.E)

                % Get indexed ellipsoid
                En = obj.E{n};

                % Get properties
                rho = En.rho;
                abc = obj.fov(:)'/2 .* En.r;
                d = obj.fov(:)'/2 .* En.d;
                A = rmat3D('zyx', pi/180*En.theta(end:-1:1));

                % Determine region and add amplitude
                ROI = vecnorm((x - d)*A' ./ abc, 2, 2) <= 1;
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