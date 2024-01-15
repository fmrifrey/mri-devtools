function data_rs = regrid(data,varargin)
% function im = recon3dflex(varargin)
%
% Part of fmrifrey/mri-devtools software package by David Frey (2023)
%   git@github.com:fmrifrey/mri-devtools.git
%
% Description: Function to quickly regrid data given zoom and
%   resolution upsampling factors
%
%
% Notes:
%   - default values in help message may not be up to date - check defaults
%       structure under the function header
%
% Path dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%
% Static Intput arguments:
%   - data:
%       - data to resample
%       - Nd float/double array with cartesian data
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'resfactor'
%       - resolution (dimension) upsampling factor
%       - float/double array describing factor along each dimension
%       - passing a value > 1 will result in higher output image dimension
%       - default is ones(ndims(data),1)
%   - 'zoomfactor'
%       - field of view zoom factor
%       - float/double array describing factor along each dimension
%       - passing a value > 1 will result in a smaller fov
%       - default is ones(ndims(data),1)
%   - 'interp'
%       - interpolation method
%       - string describing interpolation method to use in
%           griddedInterpolant()
%       - default is 'cubic'
%   - 'extrap'
%       - extrapolation method
%       - string describing extrapolation method to use in
%           griddedInterpolant()
%       - default is 'none'
%
% Function output:
%   - data_rs:
%       - resampled data
%

    % Define default arguments
    defaults = struct( ...
        'zoomfactor',   ones(1,ndims(data)), ...
        'resfactor',    ones(1,ndims(data)), ...
        'interp',       'linear', ...
        'extrap',       'none');
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});
    
    % Pad undefined dimensions with 1
    if size(args.zoomfactor(:),1) == 1
        args.zoomfactor = args.zoomfactor*ones(1,ndims(data));
    elseif size(args.zoomfactor(:),1) > ndims(data)
        error('length(zoomfactor) must be <= ndim(data)');
    elseif size(args.zoomfactor,2) < ndims(data)
        args.zoomfactor = [args.zoomfactor(:)',ones(1,ndims(data)-size(args.zoomfactor,2))];
    end
    if size(args.resfactor(:),1) == 1
        args.resfactor = args.resfactor*ones(1,ndims(data));
    elseif size(args.resfactor(:),1) > ndims(data)
        error('length(resfactor) must be <= ndim(data)');
    elseif size(args.resfactor,2) < ndims(data)
        args.resfactor = [args.resfactor(:)',ones(1,ndims(data)-size(args.resfactor,2))];
    end
    
    % Make grid arrays using imgrid
    G_in = cell(1,ndims(data));
    G_out = cell(1,ndims(data));
    [G_in{:}] = imgrid(1,size(data));
    [G_out{:}] = imgrid(1./args.zoomfactor,size(data).*args.resfactor);
    
    % Make interpolant and save resampled data
    F = griddedInterpolant(G_in{:},data,args.interp,args.extrap);
    data_rs = F(G_out{:});

end

