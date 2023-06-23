function R = rmat3D(axis,angle)
% function R = rmat3D(axis,angle)
%
% Part of fmrifrey/mri-devtools software package by David Frey (2023)
%   git@github.com:fmrifrey/mri-devtools.git
%
% Description: Function to generate rotation matrices in R3 space
%
%
% Static input arguments:
%   - axis:
%       - axis (or axises in order of mat mult.) to perform rotation on
%       - char or char string with axises in order of matrix multiplication
%       - no default, required argument
%   - angle:
%       - euler angles (rad) to perform rotation
%       - numeric value or array representing euler angles (rad)
%       - no default, required argument
%
% Function output:
%   - R:
%       - rotation matrix
%       - 3x3 matrix
%

    if numel(axis) ~= numel(angle)
        error('axises and angles must be of the same length');
    end

    % Recurse for multiple angles
    if length(axis) > 1
        R = rmat3D(axis(1:end-1),angle(1:end-1))* ...
            rmat3D(axis(end),angle(end));
    else
        switch lower(axis)
            case 'x' % X rotation matrix
                R = [
                    1,              0,              0;
                    0,              cos(angle),     -sin(angle);
                    0,              sin(angle),     cos(angle)
                    ];
            case 'y' % Y rotation matrix
                R = [
                    cos(angle),     0,              sin(angle);
                    0,              1,              0;
                    -sin(angle),    0,              cos(angle)
                    ];
            case 'z' % Z rotation matrix
                R = [
                    cos(angle),     -sin(angle),    0;
                    sin(angle),     cos(angle),     0;
                    0,              0,              1
                    ];
            otherwise
                error('invalid axis %s (argument 1)', axis);
        end
    end

end

