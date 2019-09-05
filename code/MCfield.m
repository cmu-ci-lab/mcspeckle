function [u]=MCfield(sigt,albedo,box_min,box_max,l,v,is_ff_l,is_ff_v,maxItr,lambda,doCBS,smpFlg,sct_type,ampfunc,ampfunc0,lmean0)
%MCFIELD MC rendering field algorithm
%
% render a speckle field for Nv viewings and Nl lights
%
% u=MCfield(sigt,albedo,box_min,box_max,l,v,is_ff_l,is_ff_v,maxItr,lambda,doCBS,smpFlg,sct_type,ampfunc,ampfunc0,lmean0)
%
% INPUT:
%     * 'sigt' - extinction coefficient.
%     * 'albedo' - chance to absorption event, number between 0 to 1. also
%     defined as sigs/sigt, where sigs is the scattering coefficient.
%     * 'box_min' - 2D or 3D vector for lower bound on box.
%     * 'box_max' - 2D or 3D vector for upper bound on box.
%     * 'l' - illumination direction or points 3xNl (for 3D) or 2xNl (for 
%     2D). can be defined also as 1xNl vector, in this case is interperated
%     as angles and converted to 2D vectors.
%     * 'v' - viewing directions 1xNv or 2xNv or 3xNv, directions or
%     points, defined the same as 'l'.
%     * 'is_ff_l' - true for illumination in far field, and false for 
%     illumination in near field.
%     * 'is_ff_v' - true for view in far field, and false for view in near 
%     field.
%     * 'maxItr' - number of iterations to run MC algorithm.
%     * 'lambda' - the wavelength.
%     * 'doCBS' - true for activating Coherent Backscattering.
%     * 'smpFlg' - sampling method for first particle. 1 for unifrom
%     distribted sampling, and 2 for exponential distribution sampling.
%     * 'sct_type' - scattering event type. 1 for isotropic, 2 for tabulated
%     amplitude function, 3 for Henyey-Greenstein (HG) function.
%     * 'ampfunc' - scattering function parameter. for amplitude function 
%     is a constructed table with the needed parameters (see
%     measuredFarField), for HG the g parameter.
%     * 'ampfunc0' - (optinal) scattering function for first scattering
%     event.
%     * 'lmean0' - (optional) direction of first scattering event
%
% OUTPUT:
%     * 'u' - rendered field in size of |Nv|x|Nl|

%% Check validity of some of the input
narginchk(14,16);

% get the dimensions size
if((numel(box_max) ~= 2 && numel(box_max) ~= 3) || ...
        (size(box_max,2) ~= 1) || (any(size(box_max) ~= size(box_min))))
    error('Invalid box size');
end

dim = size(box_min,1);

% get number of sources
if(size(l,1) ~= dim)
    if(dim == 2 && size(l,1) == 1)
        l = [sin(l); cos(l)];
    else
        error('Invalid light source input');
    end
end

if(size(v,1) ~= dim)
    if(dim == 2 && size(v,1) == 1)
        v = [sin(v); cos(v)];
    else
        error('Invalid view source input');
    end
end

Nl = size(l,2);
Nv = size(v,2);

%% Prepare for algorithm

% Initiate output parameters
u = zeros(Nv,Nl);

% first scattering event direction
if ~exist('lmean0','var')
    unmeanl = mean(l,2);
else
    unmeanl = lmean0;
end
meanl = unmeanl/norm(unmeanl);

% Box size
box_w = box_max-box_min;

% Pre-calculate single scattering rotation amplitude, only possible when
% both light and view are far field (otherwise it also dependent on the
% first scatter position)
af_ang_vl = zeros(Nv, Nl);
if (is_ff_v && is_ff_l)
    for j=1:Nl
        af_ang_vl(:,j)=evalampfunc_general(l(:,j)'*v,sct_type,ampfunc,dim);
    end
end

% in far field, the entrance direction to the box is fixed
if is_ff_v
    rv=v;
end

if is_ff_l
    rl=l;
end

ff_sign=-2*(is_ff_v)+1;

% threshold to begin kill particles with low weight
killThr=0.2;

%% Begin the main loop
for itr=1:maxItr
   %itr
    % Sample the first scatter
    % x: first scattering point
    % px: probability by which first point was sampled. Needed so that
    % importance sampling integration is weighted properly
    switch smpFlg
        case 1
            % uniform distribution
            x=rand(dim,1).*(box_w)+box_min; px=1;
        case 2
            % exponential distribution
            [x,px]=expSmpX(box_min,box_max,unmeanl,sigt);
    end
    %x=[0;30;0.0000001];
    %x=[-30;0;50];
    %x
    % entrance directions for near-field sources
    if ~is_ff_v
        rv=v-repmat(x,1,Nv);
        rv=rv./repmat(sum(rv.^2,1).^0.5,dim,1);
    end
    if ~is_ff_l
        rl=repmat(x,1,Nl)-l;
        rl=rl./repmat(sum(rl.^2,1).^0.5,dim,1);
    end
    
    % single scattering rotation amplitude
    if ~(is_ff_v && is_ff_l)
        for j=1:Nl
            af_ang_vl(:,j)=evalampfunc_general(rl(:,j)'*rv,sct_type,ampfunc,dim);
        end
    end
    
    % First scattering direction
    % w - sampled direction.
    % w0p - probability of the sampled direction, needed to compute inportance sampling integral correctly. 
    if ~exist('ampfunc0','var')
        w=randn(dim,1); w=w/norm(w);
        w0p=1/sqrt(2^(dim-1)*pi);
    else
        w=smpampfunc_general(meanl, sct_type,ampfunc0);     
        w0p=(evalampfunc_general(meanl'*w,sct_type,ampfunc0,dim));
    end

    % rotation due to first scattering in multiple scattering case (s
    % function in article).
    af_l=evalampfunc_general((w'*rl),sct_type,ampfunc,dim)./w0p;
    
    % complex transmission (xi function in article) and attenuation term 
    % in multiple scattering (tau function in article) between first scattering
    % particle and the light source 
    e_l0=evalphaseatt(x,l,is_ff_l,sigt,lambda,box_min,box_max);
    
    % complex volumetric throughput (ni function in article) of first
    % scattering event to be used with paths of length >1 (the multiple scattering
    % case)
    e_l0_ms=e_l0.*af_l;
    
    % complex volumetric throughput connecting first scattering particle
    % and the sensors
    e_v0=evalphaseatt(x,ff_sign*v,is_ff_v,sigt,lambda,box_min,box_max);
    
    % in case of coherent backscattering, calculate also the complex
    % volumetric throughput where the path begins from the view to the
    % first scatter
    if doCBS
         af_v=evalampfunc_general((-w'*rv),sct_type,ampfunc,dim)./w0p;
         e_v0_ms=e_v0.*af_v;
    end
      
    % number of scattering events
    pL=0;
    
    % intensity loss due to albedo
    weight=albedo;
    
    % begin paths sampling loop
    while 1

        pL=pL+1;

        % calculate the complex volumetric throughput for the last
        % scattering event in case of multiple scattering
        if (pL>1)
            e_v=evalphaseatt(x,ff_sign*v,is_ff_v,sigt,lambda,box_min,box_max);
            af_v=evalampfunc_general((ow'*rv),sct_type,ampfunc,dim);
            e_v_ms=e_v.*af_v;
            if doCBS
                e_l=evalphaseatt(x,l,is_ff_l,sigt,lambda,box_min,box_max);
                af_l=evalampfunc_general((-ow'*rl),sct_type,ampfunc,dim);
                e_l_ms=e_l.*af_l;
            end
        end

        % Update field with next-event estimation
        if (pL==1)
            tpath=af_ang_vl.*(e_v0(:)*conj(e_l0(:))');
        else
            tpath=(e_v_ms(:)*conj(e_l0_ms(:))');
            if doCBS
                tpath=1/sqrt(2)*(tpath+(e_v0_ms(:)*conj(e_l_ms(:))'));
            end
        end

        % weight path
        tpath=sqrt(weight./px)*tpath;
        % sample random phase for path
        tpath=tpath*exp(2*pi*1i*rand);
        %add path to field
        u=u+tpath;

        % advance to the next scattering event
        d=-log(-rand+1)/(sigt);
        x=x+d*w;

        % move to the next particle if the next scattering event is outside
        % the box
        if(max(x>box_max) || max(x<box_min))
            break
        end

        % albedo intensity reduction. If the intensity is too low, it is
        % possible that the particle will be killed
        if weight<killThr
            if rand>albedo
                break
            end
        else
            weight=weight*albedo;
        end

        % Sample new scatteing direction
        ow=w;
        w=smpampfunc_general(ow, sct_type,ampfunc);

    end
  
end

%% Normalization
V=prod(box_w);
u=u*sqrt(1/maxItr*V*sigt);

