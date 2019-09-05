function [scatter] = HGScatter(g, g0, albedo)
%HGSCATTER define Henyey-Greenstein scatter propertites
%
%   scatter = HGScatter(g) HG scatter with parameter g. g is scalar between
%   -1 to 1. negative g for back scattering. g = 0 for isotropic
%   scattering.
%
%   scatter = HGScatter(g,g0) g0 is a parameter spasifying how to sample
%   the first direction of a path
%   Default value is isotropic scattering (g=0).
%
%   scatter = HGScatter(g,g0,albedo) define albedo parameter. scalar from 0
%   to 1. 1 is without absorption, and 0 is full absorption. default value
%   is 1.
%
%   Class support for all inputs:
%      float: double
%
% SEE ALSO: isotropicScatter, tabulatedAmplitudeScatter, scmc
%

% Check input num
narginchk(1,5);

scatter.type = 'scatter';
scatter.function = 'HG';

if(nargin >= 3)
    if(~isscalar(albedo) || albedo < 0 || albedo > 1)
        error('Invalid albedo input');
    end
    
    scatter.albedo = albedo;
else
    % default albedo value
    scatter.albedo = 1;
end

if(nargin >= 2)
    if(~isscalar(g0) || g0 < -1 || g0 > 1)
        error('Invalid g0 input');
    end
    
    scatter.g0 = g0;
else
    % default threshold value (inf for default)
    scatter.g0 = inf;
end

if(~isscalar(g) || g < -1 || g > 1)
    error('Invalid g input');
end

scatter.g = g;

end

