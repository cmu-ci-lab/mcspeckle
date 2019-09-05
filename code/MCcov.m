function [Ms,Mm]=MCcov(sigt,albedo,box_min,box_max,l,v,is_ff_l,is_ff_v,maxItr,lambda,doCBS,smpFlg,sct_type,ampfunc,ampfunc0)
%MCCOV MC rendering of covariance algorithm
%
% Calculate the speckle covaraince of far field or near field views and
% lights.
%
% [Ms,Mm]=MCcov(sigt,albedo,box_min,box_max,l,v,is_ff_l,is_ff_v,maxItr,lambda,doCBS,smpFlg,sct_type,ampfunc,ampfunc0,lrad)
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
%     distribted sampling, and 2 for exponential sampling.
%     * 'sct_type' - scattering event type. 1 for isotropic, 2 for
%     tabulated amplitude function, 3 for Henyey-Greenstein (HG) function.
%     * 'ampfunc' - scattering function parameter. for amplitude function 
%      a constructed table with the needed parameters (see
%     measuredFarField), for HG the g parameter.
%     * 'ampfunc0' - (optinal) scattering function for first scattering
%     event.
%
% OUTPUT:
%     * 'Ms' - speckle covaraince matrix due to single scattering. 4D 
%     matrix of size |Nv|x|Nv|x|Nl|x|Nl|.
%     * 'Mm' - speckle covaraince matrix due to multiple scattering. 4D 
%     matrix of size |Nv|x|Nv|x|Nl|x|Nl|. The total covariance matrix
%     is Ms + Mm.

%% Check validity of some of the input
narginchk(14,15);

% get the dimensions 
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
Ms = zeros(Nv,Nv,Nl,Nl);
Mm = zeros(Nv,Nv,Nl,Nl);

% Scattering direction of first scattering event
unmeanl = mean(l,2);
meanl = unmeanl/norm(unmeanl);

% Box size
box_w = box_max - box_min;

% Pre-calculate single scattering rotation amplitude, only possible when
% both light and view are far field (otherwise it also dependent on the
% first scatter position)
af_ang_vl = zeros(Nv, Nl);
if (is_ff_v && is_ff_l)
    for j=1:Nl
        af_ang_vl(:,j) = evalampfunc_general(l(:,j)'*v,sct_type,ampfunc,dim);
    end
end
 
% in far field, the entrance direction to the box is fixed
if is_ff_v
    rv = v;
end

if is_ff_l
    rl = l;
end
 
ff_sign=-2*(is_ff_v)+1;

% threshold to begin kill particles with low weight
killThr=0.2;
%% Begin the main loop
for itr=1:maxItr
    
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
    % scattering event to be used with  paths of length >1 ( the multiple scattering
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
      
    % all possible volumetric throughput pair of viewing directions 
    % to be used when path length=1
    % Below we will define the matrices M00 M0n Mn0 Mnn corresponding to
    % all 4 combinations of forward and reversed paths. 
    % M00: 2 reversed paths connecting the first node on the path to v1 v2
    % M0n: reversed +forward, first node connected to v1, last node to v2
    % Mn0: forward + reversed, last node connected to v1, first node to v2
    % Mnn: forward +forward, last node connected to v1  and v2
    
    M00=e_v0(:)*e_v0(:)';
  
    % number of scattering events
    pL=0;
    
    % intensity loss due to albedo
    weight=albedo;
    
    % begin paths sampling loop
    while 1

        pL=pL+1;

        % calculate the complex volumetric throughput for the last
        % scattering event in case of path length >1
        if (pL>1)

            e_v=evalphaseatt(x,ff_sign*v,is_ff_v,sigt,lambda,box_min,box_max);
            af_v=evalampfunc_general((ow'*rv),sct_type,ampfunc,dim);
            e_v_ms=e_v.*af_v;

            if doCBS

                e_l=evalphaseatt(x,l,is_ff_l,sigt,lambda,box_min,box_max);
                af_l=evalampfunc_general((-ow'*rl),sct_type,ampfunc,dim);
                e_l_ms=e_l.*af_l;

                Mn0=e_v_ms(:)*e_v0_ms(:)';
                M0n=e_v0_ms(:)*e_v_ms(:)';
            end

            Mnn=e_v_ms(:)*e_v_ms(:)';
        end
        
        
        % in cbs, the complex volumetric throughput of the views in case of
        % the reverse path can be pre-calculated
        if ((pL==2) && (doCBS))
            M00=e_v0_ms(:)*e_v0_ms(:)';
        end

        % Update covariance with next-event estimation
        if (pL==1)
            % speckle covaraince matrix due to single scattering
            for j1=1:Nl
                for j2=1:Nl
                    Ms(:,:,j1,j2)=Ms(:,:,j1,j2)+...
                       weight/px*M00 * ...
                       (conj(e_l0(j2))*e_l0(j1)) .* ...
                       (af_ang_vl(:,j1)*af_ang_vl(:,j2)');
                end
            end
        else
            % speckle covaraince matrix due to multiple scattering
            if doCBS
                for j1=1:Nl
                    for j2=1:Nl
                        Mm(:,:,j1,j2)=Mm(:,:,j1,j2)+...
                          weight/px*(...
                          M00*(conj(e_l_ms(j2))*e_l_ms(j1))+...
                          Mnn*(conj(e_l0_ms(j2))*e_l0_ms(j1))+...
                          Mn0*(conj(e_l_ms(j2))*e_l0_ms(j1))+...
                          M0n*(conj(e_l0_ms(j2))*e_l_ms(j1)));
                    end
                end
            else
                for j1=1:Nl
                    for j2=1:Nl
                        Mm(:,:,j1,j2)=Mm(:,:,j1,j2)+...
                          weight/px*(Mnn*(conj(e_l0_ms(j2))*e_l0_ms(j1)));

                    end
                end
            end
        end

        % advance to the next scattering event
        d=-log(-rand+1)/(sigt);
        x=x+d*w;

        % Terminate the path if the next scattering event that was sampled
        % is outside the box
        if(max(x>box_max) || max(x<box_min))
            break
        end

        % albedo intensity reduction. if the intensity is too low, it is
        % possible that the particle will be killed
        if weight<killThr
            if rand>albedo
                break
            end
        else
            weight=weight*albedo;
        end
        
        % sample new scatteing direction
        ow=w;
        w=smpampfunc_general(ow, sct_type,ampfunc);
    end
end

%% Normalization
if doCBS
   Mm=Mm/2;
end

V=prod(box_w);
Mm=Mm/(maxItr)*V*sigt;
Ms=Ms/(maxItr)*V*sigt;

