function [e,nx]=evalphaseattHetro(x,v,is_ff,sigt,lambda,box_min,box_max,box_bin)



Nv=size(v,2);
dim=size(v,1);





if (is_ff)
        pv=x'*v;
        rv=ones(1,Nv);
       % bdv=cubeDist(x,box_min,box_max,-v);
       nd=10^12*ones(1,Nv);
else
    d=repmat(x,1,Nv)-v;
    nd=sum(d.^2,1).^0.5;
    v=d./repmat(nd,dim,1);
    pv=nd; 
    rv=1./(nd+0.01*lambda).^((dim-1)/2);
    %bdv=min(cubeDist(x,box_min,box_max,-v),nd);  
end
 

bdv=zeros(1,Nv);
wbdv=zeros(1,Nv);
for j=1:Nv
    
    nx=x;
    
    while 1
        [bin_ind,bin_min,bin_max]=findBin(nx,box_min,box_max,box_bin);
        [d,nx]=cubeDist(nx,bin_min,bin_max,-v(:,j));
        bdv(j)=bdv(j)+d;
        wbdv(j)=wbdv(j)+d*sigt(bin_ind);
        if (bdv(j)>=nd(j))
            md=bdv(j)-nd(j);
            bdv(j)=nd(j);
            wbdv(j)=wbdv(j)-md*sigt(bin_ind);
            break
        end
        
        if(max(nx>box_max)|max(nx<box_min))
            
            break
        end
    end
    
end
e=exp(i*pi*2/lambda*pv-1/2*wbdv).*rv;