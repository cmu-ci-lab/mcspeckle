function [u]=MCfield(sigt,albedo,box_min,box_max,l,v_max,v_stp,maxItr,lambda,doCBS,smpFlg,lmean0,vsign)
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
%
%    v_max, v_stp: (scalars) define a grid of viewing directions
%    vsign: (-1 or 1) whatever viewing grid at the front of back side of the phantom
%     
%   
%    
%     * 'maxItr' - number of iterations to run MC algorithm.
%     * 'lambda' - the wavelength.
%     * 'doCBS' - true for activating Coherent Backscattering.
%     * 'smpFlg' - sampling method for first particle. 1 for unifrom
%     distribted sampling, and 2 for exponential distribution sampling.
%     * 'sct_type' - this code is only for the isotropic scattering case
%     * 'lmean0' - (optional) direction of first scattering event



%
% OUTPUT:
%     * 'u' - rendered field in size of |Nv|x|Nl|

%% Check validity of some of the input


% get the dimensions size
if((numel(box_max) ~= 2 && numel(box_max) ~= 3) || ...
        (size(box_max,2) ~= 1) || (any(size(box_max) ~= size(box_min))))
    error('Invalid box size');
end

dim = size(box_min,1);



Nl = size(l,2);
%Nv = size(v,2);

%% Prepare for algorithm

% Initiate output parameters


% first scattering event direction
if ~exist('lmean0','var')|isempty(lmean0)
    unmeanl = [0;0;1];
else
    unmeanl = lmean0;
end
meanl = unmeanl/norm(unmeanl);

if ~exist('vsign','var')  %this input parameter indicates from which side of the phantom the viewing grid is located
    vsign=-1;
end
% Box size
box_w = box_max-box_min;






ff_sign=-2*(1)+1;


vbin_scl=1;
vrng_scl=1%5;
theta_v0=[-v_max:v_stp/vbin_scl:v_max];
v_max=v_max*vrng_scl;
theta_v=[-v_max:v_stp/vbin_scl:v_max];



Nv=length(theta_v)
Nv0=length(theta_v0)


box_stp=lambda/(v_max*2+v_stp);

Nx=ceil(box_w(1)/box_stp/2)*2+1;
box_w_n=Nx*box_stp;
box_min(1:2)=box_min(1:2)-(box_w_n-box_w(1))/2;
box_max(1:2)=box_max(1:2)+(box_w_n-box_w(1))/2;
box_w(1:2)=box_w_n;

[gx,gy]=ndgrid([box_min(1)+box_stp/2:box_stp:box_max(1)],[box_min(2)+box_stp/2:box_stp:box_max(2)]);


Nx_padd_size=max((Nv-Nx)/2,0);
Nv_padd_size=max((Nx-Nv)/2,0);
Nv0_padd_size=max((Nv-Nv0)/2,0);
Nva_padd_size=Nv0_padd_size+Nv_padd_size;

[vx,vy]=ndgrid(-theta_v,-theta_v);
vz=vsign*sqrt(1-vx.^2-vy.^2);

vx0=vx(1+Nv0_padd_size:end-Nv0_padd_size,1+Nv0_padd_size:end-Nv0_padd_size);
vy0=vy(1+Nv0_padd_size:end-Nv0_padd_size,1+Nv0_padd_size:end-Nv0_padd_size);
vz0=vz(1+Nv0_padd_size:end-Nv0_padd_size,1+Nv0_padd_size:end-Nv0_padd_size);


z_bin=lambda/max(1-abs(vz(:)))/5;
z_bin=z_bin/2;%/5;
gz=[box_min(end):z_bin:box_max(end),box_max(end)];
%Nz=length(zL);

Nz=length(gz);


diratt0=1./abs(vz0).*(vz0<0);
idiratt0=1./abs(vz0).*(vz0>0);   

l_diratt=1./abs(l(end)).*(l(end)>0);
l_idiratt=1./abs(l(end)).*(l(end)<0);   

%zexp0=(-2*pi*i/lambda*vz0-sigt/2*diratt0);
%boxexp0=exp(-sigt/2*(-box_min(end)).*diratt0);

%zexp0=(-2*pi*i/lambda*vz0-sigt/2*diratt0);
zexp0=(-2*pi*i/lambda*vz0*vsign-sigt/2*(diratt0-idiratt0));

%boxexp0=exp(-sigt/2*(-box_min(end)).*diratt0);
boxexp0=exp(-sigt/2*(-box_min(end).*diratt0+box_max(end).*idiratt0));
boxexp0=boxexp0.*exp(2*pi*i/lambda*(-box_min(end).*(vz0<0).*vz0+box_max(end).*(vz0>0).*vz0 ));




% threshold to begin kill particles with low weight
killThr=0.2;

ExpOD=sigt*box_w(end)+1;
xL=zeros(dim,round(maxItr*ExpOD));
x1L=zeros(dim,round(maxItr*ExpOD));

pxL=zeros(1,size(xL,2));
pathL=zeros(1,size(xL,2));
c_itr=0;
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
      
        x1=x;
  
        w=randn(dim,1); w=w/norm(w);
        w0p=1/sqrt(2^(dim-1)*pi);
        
        weight=albedo;
        p_len=0;
        while 1
            p_len=p_len+1;
            c_itr=c_itr+1;
            xL(:,c_itr)=x;
            x1L(:,c_itr)=x1;
            pxL(:,c_itr)=sqrt(weight./px)*exp(2*pi*i*rand);
            if (doCBS)&(p_len>1), pxL(:,c_itr)= 1/sqrt(2)*pxL(:,c_itr); end
            pathL(:,c_itr)=p_len;
            
            
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
        w=randn(dim,1);  w=w/norm(w); %smpampfunc_general(ow, sct_type,ampfunc);
        
      
        
        end
end
maxItr=c_itr;
xL=xL(:,1:maxItr);
x1L=x1L(:,1:maxItr);
pxL=pxL(:,1:maxItr);
pathL=pathL(:,1:maxItr);
 


u=zeros(Nv0,Nv0,Nl);

for tt=1:2
    
    if (tt==2)
        if ~doCBS
            continue
        end
        txL=xL;
        xL=x1L; x1L=txL;
    end
    
    [sv,si]=sort(x1L(end,:));
    xL=xL(:,si);
    pxL=pxL(:,si);
    x1L=x1L(:,si);
    pathL=pathL(:,si);
    
    
    
    z_ind_max=0;
    plane_count=zeros(1,Nz-1);
    time_count=zeros(1,Nz-1);
    for j1=1:Nz-1
        j1/Nz
        c_itr=z_ind_max+1;
        z_max=gz(j1+1);
        while (c_itr<=maxItr)& (x1L(end,c_itr)<z_max)
            c_itr=c_itr+1;
        end
        z_ind_min=z_ind_max+1;
        z_ind_max=c_itr-1;
        
        if z_ind_max<z_ind_min
            continue
        end

        jj=[z_ind_min:z_ind_max];
        tz=mean(gz(j1:j1+1));
        plane_count(j1)=length(jj);
            
        xy_inds=round((x1L(1:2,jj)-[gx(1);gy(1)])/box_stp)+1;
        inds=xy_inds(1,:)+(xy_inds(2,:)-1)*Nx;
       
       
        %parfor j2=1:Nl
        for j2=1:Nl
            
            u_xy=zeros(Nx,Nx);
            for j3=1:length(jj)
               u_xy(inds(j3))=u_xy(inds(j3))+pxL(1,jj(j3)).*exp(-2*pi*i/lambda*(-x1L(3,jj(j3))+tz)).*exp(2*pi*i/lambda*l(:,j2)'*xL(:,jj(j3))-sigt/2*((xL(end,jj(j3))-box_min(end)).*l_diratt+(box_max(end)-xL(end,jj(j3))).*l_idiratt));
            end
            

            u_xy=padarray(u_xy,[1,1]*Nx_padd_size,'both');

            tfu=fliplr(flipud(fftshift(fft2(ifftshift((u_xy))))));

            tfu=tfu(Nva_padd_size+1:end-Nva_padd_size,Nva_padd_size+1:end-Nva_padd_size);
            tfu=tfu.*exp(zexp0*tz);
            u(:,:,j2)=u(:,:,j2)+tfu;
          
        end
        
    end
    
    
end

for j2=1:Nl
    u(:,:,j2)=u(:,:,j2).*boxexp0;
end



u=u(1:vbin_scl:end,1:vbin_scl:end,:);



%% Normalization
V=prod(box_w);
u=u*sqrt(1/maxItr*V*sigt);

