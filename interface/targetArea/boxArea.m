function [targetArea] = boxArea(wavelength, MFP, z, x, y)
%BOXAREA create target area in a shape of 2D rectangle or 3D square box.
%
%   targetArea = boxArea(wavelength,scatteringMFP,z,x) defines the
%   wavelength,Mean Free Path, and diamentions of 2D rectangle target area.
%   * z is the depth of the rectangle, defined as 2 elements vector, where
%   the rectangle depth is defined from z(1) to z(2).
%   * x is the width of the rectangle, defined as 2 elements vector, where
%   the rectangle width is defined from x(1) to x(2).
% 
%   targetArea = boxArea(wavelength,scatteringMFP,z,x,y) defines 3D square
%   box, where the height defined from y(1) to y(2).
%
%   Class support for inputs wavelength, MFP:
%      float: double scalars
%
%   Class support for inputs z, x, y:
%      float: double
%
% SEE ALSO: mccov
%

% Check input num
narginchk(4,5);

if(~isscalar(wavelength) || ~isscalar(MFP))
    error('wavelength and MFP are scalars');
end

if(numel(z) ~= 2 || numel(x) ~= 2 || (nargin == 5 && numel(y) ~= 2))
    error('x, y and z defined as 2 elements vectors');
end

if(~isreal(z) || ~isreal(x) || (nargin == 5 && ~isreal(y)))
    error('x, y and z are real values');
end

if(z(1) > z(2))
    error('z(1) must be smaller than z(2)');
end

if(x(1) > x(2))
    error('x(1) must be smaller than x(2)');
end

if(nargin == 5 && y(1) > y(2))
    error('y(1) must be smaller than y(2)');
end

targetArea.type = 'targetArea';
targetArea.wavelength = wavelength;
targetArea.MFP = MFP;
targetArea.z = z;
targetArea.x = x;

if(nargin == 4)
    targetArea.D = 2;
end

if(nargin == 5)
    targetArea.D = 3;
    targetArea.y = y;
end

end
