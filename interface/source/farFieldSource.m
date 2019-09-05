function [source] = farFieldSource(in1, ~)
%FARFIELD build far field sources
%
%   source = farFieldSource(directions) N directions defined explicitly.
%   directions is of size  [2xN] or [3xN].
%     * forward scattering applies when the signs of the z component of 
%           the illumination and viewing directions are identical 
%     * backward scattering applies when the signs of the z component of 
%           the illumination and viewing directions are different 
%   source = farFieldSource(theta,0) N directions, where theta defined as
%   [1xN] or [Nx1] in radians.
%
%   Class support for all inputs:
%      float: double
%
% SEE ALSO: nearFieldSource, scmc
%

% Check input num
narginchk(1,2);

source.type = 'source';

source.farField = 1;

if(nargin == 1)
    % source position are already given
    if(size(in1,1) ~=3 && size(in1,1) ~=2)
        error('Invalid positions elements number');
    end
    
    source.count = size(in1,2);
    source.directions = in1;
        
    % normilize
    vectorLength = sqrt(sum(abs(source.directions).^2));
    source.directions = source.directions ./ vectorLength;
end

if(nargin == 2) 
    if( ~isvector(in1) )
        error('theta must be vector of directions');
    end
    
    source.count = numel(in1);
    
    source.directions2D = [sin(in1(:)), ...
            cos(in1(:))].';
    
    source.directions3D = [sin(in1(:)), ...
            zeros(source.count,1), ...
            cos(in1(:))].';
    
end


end

