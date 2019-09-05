function scatter = tabulatedAmplitudeScatter(theta,f,albedo)
%TABULATEDFIELDSCATTER define scatter with complex amplitude function.
%
%   scatter = tabulatedFieldScatter(theta,f) Define the amplitude function
%   f in directions of theta. theta is vector in size of 1xN and f is also
%   vector in size of 1xN, where N is the number of entries in directions.
%     * 2D space: theta is defined from 0 to 2pi, where 0 is the forward
%     direction of the field, and pi is the backward direction
%     * 3D space: theta is the elevation, defined from 0 to pi. 0 is the
%     forward direction and pi is the backward. the azimuth is uniformly
%     distributed from 0 to 2pi.
%
%   scatter = tabulatedFieldScatter(theta,f,albedo) define albedo
%   parameter. scalar from 0 to 1. 1 is without absorption, and 0 is full 
%   absorption. default value is 1.
%
%   More About: 
%     * Amplitude Function: the amplitude function f is measurement of the
%       complex field due to incident wave interaction with a single
%       scatterer. Those measurements can be produced with exact wave
%       solver, such as mu-diff. The amplitude function is also the square
%       root of the probability density function (pdf) of the scattering.
%
%     * Cross Section: the cross section of scatterer is defined in 2D as 
%       mean(abs(f).^2)*2*pi, and in 3D as
%       mean((abs(f).^2).*sin(theta))*2*pi^2.
%
%     * Evaluation and Sampling: The evaluation in a specific direction is
%       normalized by the cross section, thus the integral of the evaluated
%       values over all directions is normalized to 1. For sampling, we
%       sample only the elevation, while the azimuth is sampled uniformly.
%       Thus, the elevation sampling is weighted by the corresponding  
%       perimeter.
%
%   Class support for all inputs:
%      float: double
%
% SEE ALSO: isotropicScatter, HGScatter, scmc
%

% Check input num
narginchk(2,3);

if(nargin == 3)
    scatter.albedo = albedo;
else
    scatter.albedo = 1;
end

scatter.type = 'scatter';
scatter.function = 'tabulated';

if(~isvector(theta) || ~isvector(f))
    error('theta and f must be vectors')
end

if(numel(theta) ~= numel(f))
    error('theta and f must have same elements number')
end

if(~ismonotonic(theta,1))
    error('theta must be monotonic increasing')
end

if( theta(1) == 0 && theta(end) == 2 * pi)
    scatter.D = 2;
else
    if(theta(1) == 0 && theta(end) == pi)
       scatter.D = 3;
    else
        error('theta must be from 0 to 2pi in 2D or -1 to 1 in 3D')
    end
end
    
scatter.ampfunc.theta = theta(:).';

% Prepare the complex amplitude function for 2D and 3D
if(scatter.D == 2)
    scatter.ampfunc.cs = mean(abs(f).^2)*2*pi;
    scatter.ampfunc.samplePdf = (abs(f(:).').^2) ./ sum((abs(f(:).').^2));
    scatter.ampfunc.sampleCdf = cumsum(scatter.ampfunc.samplePdf);
    scatter.ampfunc.sampleIcdf = invCDF(scatter.ampfunc.sampleCdf,2);
    
    scatter.ampfunc.evalPdf = (abs(f(:).').^2) ./ scatter.ampfunc.cs;
    scatter.ampfunc.evalAmp = scatter.ampfunc.evalPdf .^ 0.5 ...
        .* exp(1i*angle(f(:).'));
end

if(scatter.D == 3)
    scatter.ampfunc.cs = mean((abs(f).^2).*sin(theta))*2*pi^2;
    scatter.ampfunc.samplePdf = (abs(f(:).').^2) .* sin(theta) ./ ...
        sum((abs(f(:).').^2) .* sin(theta));
    scatter.ampfunc.sampleCdf = cumsum(scatter.ampfunc.samplePdf);
    scatter.ampfunc.sampleIcdf = invCDF(scatter.ampfunc.sampleCdf,3);
    
    scatter.ampfunc.evalPdf = (abs(f(:).').^2) ./ scatter.ampfunc.cs;
    scatter.ampfunc.evalAmp = scatter.ampfunc.evalPdf .^ 0.5 ...
        .* exp(1i*angle(f(:).'));
end

end
