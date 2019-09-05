function [scatter] = isotropicScatter(albedo)
%HGSCATTER define isotropic scattering for scatter
%
%   scatter = isotropicScatter() isotropic scattering.
%
%   scatter = isotropicScatter(albedo) define albedo parameter. scalar from
%   0 to 1. 1 is without absorption, and 0 is full absorption. default
%   value is 1.
%
%   Class support for all inputs:
%      float: double
%
% SEE ALSO: HGScatter, tabulatedAmplitudeScatter, scmc
%

% Check input num
narginchk(0,1);

scatter.type = 'scatter';
scatter.function = 'isotropic';

if(nargin == 1)
    if(~isscalar(albedo) || albedo < 0 || albedo > 1)
        error('Invalid albedo input');
    end
    
    scatter.albedo = albedo;
else
    % default albedo value
    scatter.albedo = 1;
end

end

