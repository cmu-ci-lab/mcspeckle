function [meanU]=evalMeanUnifCtr(box_min,box_max,l,v,sigt,lambda)

%This is a mean calculation that is valid in the far field only, for
%homogenious case, and assuming the illumination and viewing directions are
%at a small angle with the normail of the target sample.

box_w=box_max-box_min;

dim=size(box_min,1);

Nl=size(l,2);
Nv=size(v,2);
meanU=zeros(Nv,Nl);
for j=1:Nl
  w=v(1:dim-1,:)-l(1:dim-1,j);
  w_e=v(dim,:)-l(dim,j);
  jj=find(abs(w_e)<0.7);
  aw=asin(w(:,jj));
  e=1;
  for j2=1:dim-1
  e=e.* sinc(aw(j2,:).*box_w(j2));
  end
  e=e*prod(box_w(1:dim-1))*exp(-sigt/2*box_w(end));
  meanU(jj(:),j)=e(:);
end