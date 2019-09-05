function [source] = nearFieldSource(in1, in2, in3)
%NEARFIELD build near field sources
%
%   source = nearFieldSource(X,Y,N) make N sources spread uniformly between point X to Y.
%     X and Y are [2x1] or [3x1] vectors, and N is scalalr.
% 
%   source = nearFieldSource(positions) make N sources  defined 
%   explicity by positions, which is [2xN] or [3xN] vector.
%
%   source = nearFieldSource(r,theta) make N sources in distance of r from
%   [0,0,0] and direction theta [1xN] or [Nx1] in radians.
%
%   Class support for all inputs:
%      float: double
%
% SEE ALSO: farFieldSource, scmc
%

% Check input num
narginchk(1,3);

source.type = 'source';

source.farField = 0;

if(nargin == 1)
    % source position are already given
    if(size(in1,1) ~=3 && size(in1,1) ~=2)
        error('Invalid positions elements number');
    end
    
    source.count = size(in1,2);
    
    source.positions = in1;
end

if(nargin == 2)
    % r is positive scalar
    if(~(isscalar(in1) && in1 > 0))
        error('r must be positive scalar');
    end
    
    if( ~isvector(in2) )
        error('theta must be vector of directions');
    end
    
    source.count = numel(in2);

    source.positions2D = in1 * [sin(in2(:)), ...
        cos(in2(:))].';
    
    source.positions3D = in1 * [sin(in2(:)), ...
            zeros(source.count,1), ...
            cos(in2(:))].';
    
end

if(nargin == 3)
    % source position are given as interpulation between two points
    if((size(in1,1) ~=3 && size(in1,1) ~=2) || size(in1,2) ~= 1)
        error('Invalid X position');
    end
    
    if((size(in2,1) ~=3 && size(in2,1) ~=2) || size(in2,2) ~= 1)
        error('Invalid Y position');
    end
    
    if(size(in1,1) ~= size(in2,1))
        error('Both points must have the same diamentions');
    end
    
    if(~(isscalar(in3) && rem(in3,1) == 0 && in3 > 0))
        error('Number of input points must be positive integer');
    end
    
    source.count = in3;
    
    if(size(in1,1) == 2)
        source.positions = [linspace(in1(1),in2(1),in3).', ...
            linspace(in1(2),in2(2),in3).'].';
    else
        source.positions = [linspace(in1(1),in2(1),in3).', ...
            linspace(in1(2),in2(2),in3).', ...
            linspace(in1(3),in2(3),in3).'].';
    end
    

end


end
