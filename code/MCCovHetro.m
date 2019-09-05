function [Ms,Mm,mean1]=MCCovHetro( sigt, albedo, box_min,box_max,box_bin, l, v,is_ff_l,is_ff_v,maxItr,lambda,doCBS,smpFlg,sct_type,ampfunc,ampfunc0,lrad)
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
%doCBS=1;
Nl=size(l,2);
Nv=size(v,2);

dim=size(box_min,1);


if size(l,1)==1
    l=[sin(l); cos(l)];
end
if size(v,1)==1
    v=[sin(v); cos(v)];
end

if (length(albedo)==1)
    albedo=albedo*ones(size(sigt));
end


Ms=zeros(Nv,Nv,Nl,Nl);
Mm=zeros(Nv,Nv,Nl,Nl);

mean1=zeros(Nv,Nl);


unmeanl=mean(l,2);
meanl=unmeanl/norm(unmeanl);
box_w=box_max-box_min;
 

if exist('lrad','var')
    lstart=(box_min+box_max)/2;
    lstart(end)=(l(end)>0)*box_min(end)+(l(end)<=0)*box_max(end);
end
 if (is_ff_v*is_ff_l)
    for j=1:Nl
      af_ang_vl(:,j)=evalampfunc_general(l(:,j)'*v,sct_type,ampfunc,dim);
    end
 end
 
 if is_ff_v
     rv=v;
 end
 if is_ff_l
     rl=l;
 end
 
 ff_sign=-2*(is_ff_v)+1;
 
 bin_prob=sigt(:)'/sum(sigt(:));
 inv_bin_prob=invCDFvec(cumsum(bin_prob),length(sigt(:))*100);
 V=sum(sigt(:))*prod(box_bin);%/mean(sigt(:));

 
 if (dim==2)
     %[gz,gx]=ndgrid([box_min(2):box_bin(2):box_max(2)],[box_min(1):box_bin(1):box_max(1)] );
     [gx,gz]=ndgrid([box_min(1):box_bin(1):box_max(1)],[box_min(2):box_bin(2):box_max(2)] );
     
     gx_min=gx(1:end-1,1:end-1);
     gz_min=gz(1:end-1,1:end-1);
     
     bin_min=[gx_min(:)';gz_min(:)'];
     bin_max=bin_min+box_bin;
 else
     [gx,gy,gz]=ndgrid([box_min(1):box_bin(1):box_max(1)],[box_min(2):box_bin(2):box_max(2)], [box_min(3):box_bin(3):box_max(3)] );
     gx_min=gx(1:end-1,1:end-1,1:end-1);
     gy_min=gy(1:end-1,1:end-1,1:end-1);
     gz_min=gz(1:end-1,1:end-1,1:end-1);
     
     bin_min=[gx_min(:)';gy_min(:)';gz_min(:)'];
     bin_max=bin_min+box_bin;
 end
 
 %x1L=zeros(dim,maxItr);
for itr=1:maxItr
    %itr/maxItr
   
    switch smpFlg
        case 1
            
            bin_ind=smpicdf(inv_bin_prob);
            
            x=rand(dim,1).*(box_bin)+bin_min(:,bin_ind); px=1/V;%sigt(bin_ind);
            if exist('lrad','var')
                r=norm(x-(l'*x)*l-lstart);
                while r>lrad
                    x=rand(dim,1).*(box_w)+box_min; px=1;
                    r=norm(x-(l'*x)*l-lstart);
                end
            end
        case 2 
            %error('not implemented yet')
            %[x,px]=expSmpX(box_min,box_max,unmeanl,sigt);
            if ~exist('lrad','var')
                lrad=[];
            end
            [x,px]=expSmpXHetro(box_min,box_max,box_bin,meanl,sigt,lrad);
            
            
    end
    if(max(x>box_max)|max(x<box_min))
        continue
    end
    
    
    %px= px*sigt(findBin(x,box_min,box_max,box_bin));
    %x1L(:,itr)=x;
    
    if ~is_ff_v
      rv=v-repmat(x,1,Nv);
      rv=rv./repmat(sum(rv.^2,1).^0.5,dim,1);
    end
    if ~is_ff_l
      rl=repmat(x,1,Nl)-l;
      rl=rl./repmat(sum(rl.^2,1).^0.5,dim,1);
    end
    if ~(is_ff_v*is_ff_l)
        for j=1:Nl
            af_ang_vl(:,j)=evalampfunc_general(rl(:,j)'*rv,sct_type,ampfunc,dim);
        end
    end
    if ~exist('ampfunc0','var')
      w=randn(dim,1); w=w/norm(w);
      w0p=1/sqrt(2^(dim-1)*pi);
    else
        
        w=smpampfunc_general(meanl, sct_type,ampfunc0);     
        w0p=(evalampfunc_general(meanl'*w,sct_type,ampfunc0,dim));
    end
  
  
         
    
    
    af_l=evalampfunc_general((w'*rl),sct_type,ampfunc,dim)./w0p;
    e_l0=evalphaseattHetro(x,l,is_ff_l,sigt,lambda,box_min,box_max,box_bin);
    e_l0_ms=e_l0.*af_l;
  
  
    u_l0=evalphaseatt(x,l,is_ff_l,0,lambda,box_min,box_max);
    e_v0=evalphaseattHetro(x,ff_sign*v,is_ff_v,sigt,lambda,box_min,box_max,box_bin);
    if doCBS
         af_v=evalampfunc_general((-w'*rv),sct_type,ampfunc,dim)./w0p;
          e_v0_ms=e_v0.*af_v;
    end
      
  if (max(isnan(e_v0))|max(isnan(e_l0)))
      e_l0
      keyboard
  end
    
    M00=e_v0(:)*e_v0(:)';
  
    
    pL=0;
    weight=albedo(findBin(x,box_min,box_max,box_bin));
  while 1
    
     pL=pL+1;
      
     if (pL>1)
         
         e_v=evalphaseattHetro(x,ff_sign*v,is_ff_v,sigt,lambda,box_min,box_max,box_bin);
         af_v=evalampfunc_general((ow'*rv),sct_type,ampfunc,dim);
         e_v_ms=e_v.*af_v;
         
         if doCBS
             
             e_l=evalphaseattHetro(x,l,is_ff_l,sigt,lambda,box_min,box_max,box_bin);
             af_l=evalampfunc_general((-ow'*rl),sct_type,ampfunc,dim);
             e_l_ms=e_l.*af_l;
             
             Mn0=e_v_ms(:)*e_v0_ms(:)';
             M0n=e_v0_ms(:)*e_v_ms(:)';
         end
         Mnn=e_v_ms(:)*e_v_ms(:)';
     end
     if ((pL==2)&(doCBS))
       M00=e_v0_ms(:)*e_v0_ms(:)';
      
     end
     
  
      if (pL==1)
                  
          mean1=mean1+weight./px*af_ang_vl.*(e_v0(:)*conj(u_l0(:))');
  
         for j1=1:Nl
             for j2=1:Nl
               Ms(:,:,j1,j2)=Ms(:,:,j1,j2)+...
                   weight/px*M00*(conj(e_l0(j2))*e_l0(j1)).*(af_ang_vl(:,j1)*af_ang_vl(:,j2)');
                if max(max(isnan(Ms(:,:,j1,j2))))
                    keyboard
                end
             end
         end
      else
          
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
      if 0
          break
      end
      %d=-log(-rand+1)/(sigt);
      d=smpWoodcock(x,box_min,box_max,box_bin,sigt,w);
      
      x=x+d*w;
      
      if(max(x>box_max)|max(x<box_min))
          break
      end
      killThr=0.2;
      if weight<killThr
          if rand>albedo(findBin(x,box_min,box_max,box_bin))
              break
          end
      else
          weight=weight*albedo(findBin(x,box_min,box_max,box_bin));
      end
      ow=w;
      w= smpampfunc_general(ow, sct_type,ampfunc);
     
  end
  
end



if doCBS
   Mm=Mm/2;
end
%V=prod(box_w);


Mm=Mm/(maxItr);
Ms=Ms/(maxItr);

if sct_type==2
    cs=mean(ampfunc.pdf)*2*pi;
else
    cs=1;
end
mean1=mean1/maxItr/sqrt(cs);
