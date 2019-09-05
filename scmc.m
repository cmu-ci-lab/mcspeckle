function [mulRes] = scmc(targetArea, views, lights, scatter, iterations, varargin)
%SCMC Speckle Covariance Monte-Catlo solver
%
%   mulRes = scmc(targetArea,views,lights,scatter,iterations) render
%   covariance of given parameters. Return structre with some of
%   the possible measurements: configuration struct, covariance,
%   single scattering covariance, rendered speckele field.
%     * 'targetArea' - define the sample properties and the wavelength.
%     see supporting functions: boxArea.
%     * 'views' - define far field views directions or near field views
%                 positions.
%     see supporting functions: farFieldSource, nearFieldSource.
%     * 'lights' - define far field lights directions or near field lights
%                  positions.
%     see supporting functions: farFieldSource, nearFieldSource.
%     * 'scatter' - define the scattering function of each scatterer. The
%                   scattering function can be isotropic, Henyey-Greenstein 
%                   (HG),and user defined tabulated function.
%     see supporting functions: isotropicScatter, HGScatter,
%                               tabulatedAmplitudeScatter.
%     * 'iterations' - number of mc iterations.
%
%   mulRes = scmc(targetArea,views,lights,scatter,iterations,varargin) add
%   options to scmc, in pair of field name and value, with the following
%   options:
%     * 'parforIters' - number of iterations using parfor. Thus the total
%                       iteartions is parforIters * iterations. The default
%                       is not using parfor.
%     * 'rng' - rng number for random number generator (positive scalar). Not
%               possible when using parfor. Default is not using rng.
%     * 'CBS' - true= evaluat coherent back scattering. Default: true.
%     * 'uniformFirstScatter' - true if the first scatterer is sampled
%                               uniformly in the target area, otherwise it
%                               is sampled according to an exponentially
%                               decaying function of the distance from the
%                               edge of the target. Default: false.
%     * 'renderCov' - true for rendering covariance.
%                     Default: true.
%     * 'singleScattering' - true for returning single scattering
%                            covariance. Default: false.
%     * 'multipleScattering' - true for returning multiple scattering
%                              covariance. Default: true.
%     * 'renderField' - true for sampling a field.
%                       Default: false.
%
%   mulRes = scmc(Config) run the algorithm with pre-calculated Config
%   struct, which is returned from scmc.m. 
%
%   Class support for targetArea, views, lights, scatter:
%      struct
%
%   Class photonsNum:
%      float: double
%
% SEE ALSO: boxArea, nearFieldSource, farFieldSource, tabulatedAmplitudeScatter, HGScatter
%

%% Check valid input
if(nargin ~= 1)
    narginchk(5,inf);

    if(mod(length(varargin),2) ~= 0)
        error('Parameters must be in struct of option and value')
    end

    if(~strcmp(targetArea.type,'targetArea'))
        error('Invalid target area input')
    end

    if(~strcmp(views.type,'source'))
        error('Invalid views input')
    end

    if(~strcmp(lights.type,'source'))
        error('Invalid lights input')
    end

    if(~strcmp(scatter.type,'scatter'))
        error('Invalid scatter input')
    end

    if(~isscalar(iterations) || iterations < 1)
        error('Invalid photons number')
    end

    % measuredFarField must be in 2D
    if(strcmp(scatter.function,'measuredFarField') && (targetArea.D == 3))
        error('Measured scattering of far field is used only in 2D')
    end
    
    % tabulated scattering need to fit to its dimensions
%     if(strcmp(scatter.function,'tabulated') && (targetArea.D ~= scatter.D))
%         error('Tabulated scattering function is not fit to the area dimensions')
%     end
    
    % scattering 2D / 3D ...
    if(targetArea.D == 3)
        if isfield(views,'directions3D')
            views.directions = views.directions3D;
        end
        if isfield(views,'positions3D')
            views.positions = views.positions3D;
        end
        
        if isfield(lights,'directions3D')
            lights.directions = lights.directions3D;
        end
        if isfield(lights,'positions3D')
            lights.positions = lights.positions3D;
        end
    end
    
    if(targetArea.D == 2)
        if isfield(views,'directions2D')
            views.directions = views.directions2D;
        end
        if isfield(views,'positions2D')
            views.positions = views.positions2D;
        end
        
        if isfield(lights,'directions2D')
            lights.directions = lights.directions2D;
        end
        if isfield(lights,'positions2D')
            lights.positions = lights.positions2D;
        end
    end
        
    %% Build Config
    Config.targetArea = targetArea;
    Config.views = views;
    Config.lights = lights;
    Config.scatter = scatter;
    Config.render.iterations = iterations;

    %% Default values
    Config.parforIters = 1;
    Config.rng = -1;
    Config.CBS = true;
    Config.uniformFirstScatter = false;
    Config.render.singleScattering = false;
    Config.render.cov = true;
    Config.render.singleScattering = false;
    Config.render.multipleScattering = true;
    Config.render.field = false;

    %% Set preffered values
    for optNum = 1:2:length(varargin)
        optString = char(varargin{optNum});
        optVal = varargin{optNum + 1};

        if strcmp(optString,'parforIters')
           if(~isscalar(optVal) || optVal < 0)
            error('invalid parfor iterations number')
           end
           
           if(Config.rng ~= -1)
               error('rng can not be used while using parfor')
           end

           Config.parforIters = optVal;
           continue;
        end
        
        if strcmp(optString,'rng')
           if(~isscalar(optVal) || optVal < 0)
            error('invalid rng number')
           end
           
           if(Config.parforIters ~= 1)
               error('rng can not be used while using parfor')
           end

           Config.rng = optVal;
           continue;
        end

        if strcmp(optString,'CBS')
           if(~isscalar(optVal) || ~islogical(optVal))
            error('CBS must be logical scalar')
           end

           Config.CBS = optVal;
           continue;
        end

        if strcmp(optString,'uniformFirstScatter')
           if(~isscalar(optVal) || ~islogical(optVal))
            error('uniformFirstScatter must be logical scalar')
           end

           Config.uniformFirstScatter = optVal;
           continue;
        end

        if strcmp(optString,'renderCov')
           if(~isscalar(optVal) || ~islogical(optVal))
            error('renderCov must be logical scalar')
           end

           Config.render.cov = optVal;
           continue;
        end

        if strcmp(optString,'singleScattering')
           if(~isscalar(optVal) || ~islogical(optVal))
            error('Single scattering option must be logical scalar')
           end

           Config.render.singleScattering = optVal;

           continue;
        end

        if strcmp(optString,'multipleScattering')
           if(~isscalar(optVal) || ~islogical(optVal))
            error('Multiple scattering option must be logical scalar')
           end

           Config.render.multipleScattering = optVal;

           continue;
        end 

        if strcmp(optString,'renderField')
           if(~isscalar(optVal) || ~islogical(optVal))
            error('renderField must be logical scalar')
           end

           Config.render.field = optVal;
           continue;
        end

        error('Invalid option value')
    end
else
    % in pre-allocated config, the first variable is the config
    Config = targetArea;
end

%% run Monte-Carlo
mulRes.Config = Config;

% build the parameters for running the algortihms
sigt = 1./Config.targetArea.MFP;

albedo = Config.scatter.albedo;

if(strcmp(Config.scatter.function,'isotropic'))
    sct_type = 1;
    ampfunc = 0;
    
    covRendParams = 14;
    fieldRenderParams = 14;
end

if(strcmp(Config.scatter.function,'tabulated'))
    sct_type = 2;
    ampfunc = scatter.ampfunc;
    
    if(isinf(Config.scatter.ampfunc0))
        covRendParams = 14;
        fieldRenderParams = 14;
    else
        covRendParams = 15;
        fieldRenderParams = 15;
        
        ampfunc0 = Config.scatter.ampfunc0;
    end
    
    
end

if(strcmp(Config.scatter.function,'HG'))
    sct_type = 3;
    ampfunc = Config.scatter.g;
    
    if(isinf(Config.scatter.g0))
        % default value to g0
        covRendParams = 14;
        fieldRenderParams = 14;
    else
        covRendParams = 15;
        fieldRenderParams = 15;
        ampfunc0 = Config.scatter.g0;
    end
end

if(Config.targetArea.D == 3)
    box_min = [Config.targetArea.x(1);Config.targetArea.y(1);Config.targetArea.z(1)];
    box_max = [Config.targetArea.x(2);Config.targetArea.y(2);Config.targetArea.z(2)];
else
    box_min = [Config.targetArea.x(1);Config.targetArea.z(1)];
    box_max = [Config.targetArea.x(2);Config.targetArea.z(2)];
end

if(Config.lights.farField == 1)
    is_ff_l = 1;
    l = Config.lights.directions;
else
    is_ff_l = 0;
    l = Config.lights.positions;
end

if(Config.views.farField == 1)
    is_ff_v = 1;
    v = Config.views.directions;
else
    is_ff_v = 0;
    v = Config.views.positions;
end

maxItr = Config.render.iterations;

lambda = Config.targetArea.wavelength;

doCBS = Config.CBS;

if(Config.uniformFirstScatter)
    smpFlg = 1;
else
    smpFlg = 2;
end

% run cov rendering
if(Config.render.cov)
    if(Config.parforIters == 1)
        if(Config.rng ~= -1)
            rng(Config.rng);
        end
        
        if(covRendParams == 14)
            [Ms,Mm] = MCcov( sigt, albedo, box_min, box_max,  ...
                l, v, is_ff_l, is_ff_v, maxItr, lambda, doCBS, smpFlg, ...
                sct_type, ampfunc);
        end

        if(covRendParams == 15)
            [Ms,Mm] = MCcov( sigt, albedo, box_min, box_max,  ...
                l, v, is_ff_l, is_ff_v, maxItr, lambda, doCBS, smpFlg, ...
                sct_type, ampfunc,ampfunc0);
        end
    else
        
        Ms = zeros(Config.views.count, Config.views.count, ...
            Config.lights.count, Config.lights.count, Config.parforIters);
        Mm = Ms;
    
        if(covRendParams == 14)
            parfor iterNum = 1:1:Config.parforIters
                [Ms(:,:,:,:,iterNum),Mm(:,:,:,:,iterNum)] = ...
                    MCcov( sigt, albedo, box_min, box_max,  l, v, ...
                    is_ff_l, is_ff_v, maxItr, lambda, doCBS, smpFlg, ...
                    sct_type, ampfunc);
            end
        end
        if(covRendParams == 15)
            parfor iterNum = 1:1:Config.parforIters
                [Ms(:,:,:,:,iterNum),Mm(:,:,:,:,iterNum)] = ...
                    MCcov( sigt, albedo, box_min, box_max,  l, v, ...
                    is_ff_l, is_ff_v, maxItr, lambda, doCBS, smpFlg, ...
                    sct_type, ampfunc, ampfunc0);
            end
        end
        
        Ms = mean(Ms,5);
        Mm = mean(Mm,5);
    end
    
    if(Config.render.singleScattering)
        mulRes.Csingle = Ms;
    end
    
    if(Config.render.multipleScattering)
        mulRes.C = Ms + Mm;
    end
    
end

% run field rendering
if(Config.render.field)
    if(Config.parforIters == 1)
        if(Config.rng ~= -1)
            rng(Config.rng);
        end
        
        if(fieldRenderParams == 14)
            u = MCfield( sigt, albedo, box_min, box_max, l, v, ...
                is_ff_l, is_ff_v, maxItr, lambda, doCBS,smpFlg, sct_type, ...
                ampfunc);
        end

        if(fieldRenderParams == 15)
            u = MCfield( sigt, albedo, box_min, box_max, l, v, ...
                is_ff_l, is_ff_v, maxItr, lambda, doCBS,smpFlg, sct_type, ...
                ampfunc, ampfunc0);
        end
    else
        u = zeros(Config.views.count, Config.lights.count, ...
            Config.parforIters);
    
        if(fieldRenderParams == 14)
            parfor iterNum = 1:1:Config.parforIters
                u(:,:,iterNum) = MCfield( sigt, ...
                    albedo, box_min, box_max, l, v, is_ff_l, is_ff_v,  ...
                    maxItr, lambda, doCBS, smpFlg, sct_type, ampfunc);
            end
        end

        if(fieldRenderParams == 15)
            parfor iterNum = 1:1:Config.parforIters
                u(:,:,iterNum) = MCfield( sigt, ...
                    albedo, box_min, box_max, l, v, is_ff_l, is_ff_v,  ...
                    maxItr, lambda, doCBS, smpFlg, sct_type, ampfunc, ampfunc0);
            end
        end
        
        % normalize
        u = u * sqrt(maxItr);
        u = sum(u,3);
        u = u./sqrt(maxItr * Config.parforIters);
        
    end
    
    mulRes.field = u;
    
end

end


