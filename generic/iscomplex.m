function TF = iscomplex(val)
% function TF = iscomplex(val)
%
% Part of fmrifrey/mri-devtools software package by David Frey (2023)
%   git@github.com:fmrifrey/mri-devtools.git
%
% Description: Quick macro to determine if value is complex
%
%
% Static input arguments:
%   - val:
%       - value (or array of values) to determine if it is complex
%       - double/float scalar or array
%       - no default, required argument
%
% Function output:
%   - TF:
%       - logical value containing answer
%       - logical scalar or array of same size as val
%

    % Initialize TF
    TF = false(size(val));
    
    % Determine if any value in val has an imaginary value
    for i = 1:length(val(:))
        TF(i) = (imag(val(i)) > 0);
    end
    
end