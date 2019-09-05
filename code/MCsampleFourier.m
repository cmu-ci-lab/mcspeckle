function [u]=MCsample( sigt, albedo, box_min,box_max, l, v_max,v_stp,maxItr,lambda,smpFlg,sct_type,ampfunc,ampfunc0,lmean0)
% is_ff_l: (binary) illumination in far field. 
% is_ff_v: (binary) viewing in far field 
% l: illumination direction 3xNl or 2xNl or 1xNl(this case is interperated
%                                               as angles and converted to
%                                               2D vectors)
%    if in far field this are assumed to be unit norm vector. In near field
%    sensor points
% v: viewing directions 1xNv or 2xNv or 3xNv, directions or points
% sct_type=1 isotropic sct_type=2 parametric amplitude function provided
% sct_type=3 HG, amp_func= g param 
% box_min, box_max: 2D or 3D vectors lower and upper bounds on box 

Nl=size(l,2);


dim=size(box_min,1);

if size(l,1)==1
    l=[sin(l); cos(l)];
end

fv_max=1/v_stp*lambda;
fv_stp=1/v_max/4*lambda;
fv_max=box_max(1);
fv_grid=[0:fv_stp:fv_max];
fv_grid=[-fv_grid(end:-1:2),fv_grid];
[fv_x,fv_y]=meshgrid(fv_grid);



 
Nv0=length(fv_grid);
Nv=Nv0^2; 
fu=zeros( Nv0,Nv0,Nl);
v0=[0;0;1]; 
%u=zeros(Nv,Nl);

if ~exist('lmean0','var')
unmeanl=mean(l,2);
else
    %lmean0
    unmeanl=lmean0;
end
meanl=unmeanl/norm(unmeanl);
box_w=box_max-box_min;
 


 for j=1:Nl
   af_ang_vl(1,j)=evalampfunc_general(l(:,j)'*v0,sct_type,ampfunc,dim);
 end
 
 
 
for itr=1:maxItr
    %itr/maxItr
   
    switch smpFlg
        case 1
            x=rand(dim,1).*(box_w)+box_min; px=1;
            
        case 2    
            [x,px]=expSmpX(box_min,box_max,unmeanl,sigt);
            
    end
    %x(1:2)=0;
    %x1L(:,itr)=x;
    
   
    if ~exist('ampfunc0','var')
      w=randn(dim,1); w=w/norm(w);
      w0p=1/sqrt(2^(dim-1)*pi);
    else
        
        w=smpampfunc_general(meanl, sct_type,ampfunc0);     
        w0p=(evalampfunc_general(meanl'*w,sct_type,ampfunc0,dim));
    end
  
  
         
    
    
    af_l=evalampfunc_general((w'*l),sct_type,ampfunc,dim)./w0p;
    e_l0=evalphaseatt(x,l,1,sigt,lambda,box_min,box_max);
    e_l0_ms=e_l0.*af_l;
  
  
    
    e_v0=evalphaseatt(x,-v0,1,sigt,lambda,box_min,box_max);
    
      
  
    
  
    
    pL=0;
    weight=albedo;
  while 1
    
     pL=pL+1;
      
     if (pL>1)
         
         e_v=evalphaseatt(x,-v0,1,sigt,lambda,box_min,box_max);
         af_v=evalampfunc_general((ow'*v0),sct_type,ampfunc,dim);
         e_v_ms=e_v.*af_v;
         
        
     end
      if (pL==1)
          
          tpath=af_ang_vl.*(e_v0(:)*conj(e_l0(:))');
          
          
      else
          tpath=(e_v_ms(:)*conj(e_l0_ms(:))');
          
        
          
      end
      tpath=sqrt(weight./px)*tpath;
      tpath=tpath*exp(2*pi*i*rand);
      %u=u+tpath;
  
      ind1=round((x(1)+fv_max)/fv_stp)+1;
      ind2=round((x(2)+fv_max)/fv_stp)+1;
      if (ind1>0)&(ind1<=Nv0)&(ind2>0)&(ind2<=Nv0)
          fu(ind2,ind1,:)=fu(ind2,ind1,:)+reshape(tpath,1,1,Nl);
      end
      
 
      
      
      d=-log(-rand+1)/(sigt);
      x=x+d*w;
      
      if(max(x>box_max)|max(x<box_min))
          break
      end
      killThr=0.2;
      if weight<killThr
          if rand>albedo
              break
          end
      else
          weight=weight*albedo;
      end
      ow=w;
      w= smpampfunc_general(ow, sct_type,ampfunc);
     
  end
  
end

for j=1:Nl
u(:,:,j)=fft2(fu(:,:,j));
end
V=prod(box_w);
u=u*sqrt(1/maxItr*V*sigt);



 