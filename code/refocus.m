function rU=refocus(u,v,l,box_min,box_max,box_stp,lambda)
%l is assumed to be a 3x1 ellements vector, v a 3xN  (that is only single
%illumination direction)
dim=size(box_min,1);

if (dim==2)
     [gx,gz]=ndgrid([box_min(1):box_stp(1):box_max(1)],[box_min(2):box_stp(2):box_max(2)] );     
     
     p_grid=[gx(:)';gz(:)'];
     
     dimv=size(gx);
     
 else
     [gx,gy,gz]=ndgrid([box_min(1):box_stp(1):box_max(1)],[box_min(2):box_stp(2):box_max(2)], [box_min(3):box_stp(3):box_max(3)] );
     
     p_grid=[gx(:)';gy(:)';gz(:)'];
    dimv=size(gx);
end
 
Ns=50;

N=size(p_grid,2);

rU=zeros(1,N);

for j=1:Ns:N
    %j/N
    jj=[j+1:min(j+Ns,N)];
    e=exp(-2*pi*i/lambda*(v'*p_grid(:,jj)));
    el=exp(2*pi*i/lambda*(l'*p_grid(:,jj)));
    rU(jj)=(u(:)'*e).*el;
end

rU=reshape(rU,dimv);
