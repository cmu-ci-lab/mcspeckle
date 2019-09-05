
lambda=1;
box_min=[-50; -50;  0];
box_max=[ 50;  50;  100];
sigt=[1/25,1/25,0,0;0,0,1/25,1/25]';
%sigt=[0,1/2000,1/2000,0;0,1/2000,1/2000,0]';
sigt=repmat(reshape(sigt,4,1,2),[1,4,1]);
box_bin=[25;25;50];

sigt=[1/50,1/25,1/25,0,0;0,0,1/50,1/25,0]';
sigt=repmat(reshape(sigt,5,1,2),[1,5,1]);
sigt(:,[1:2,4:5],:)=0;
box_bin=[20;20;50];


theta_l=pi-deg2rad([-35:5:35]);
theta_v=deg2rad([-40:0.2:40]);


Nl0=length(theta_l);
Nv0=length(theta_v);
Nv=Nv0^2;
Nl=Nl0^2;

Nx=(box_max(1)-box_min(1))/box_bin(1);
Nz=(box_max(2)-box_min(2))/box_bin(2);
Nxz=Nx*Nz;

[vx,vy]=ndgrid(theta_v,theta_v);
v=[sin(vx(:))';sin(vy(:))';ones(1,Nv)];
v=v./repmat(sum(v.^2,1).^0.5,3,1);
[lx,ly]=ndgrid(theta_l,theta_l);
l=[sin(lx(:))';sin(ly(:))';ones(1,Nl)];
l=l./repmat(sum(l.^2,1).^0.5,3,1);




doCBS=1; smpFlg=1;
maxItr=10^4;
%maxItr=1;
maxItr=200;
j=1
rng(520)
parfor j=1:24
tic
uL(:,:,j)=MCsampleHetro( sigt, 1, box_min,box_max,box_bin, l, v,1,1,maxItr,lambda,doCBS,smpFlg,1,0);
toc
end
uL=sum(uL,3);
box_stp=[1;50;2];
box_stp=[2;50;5];

j=ceil(Nl/2);
trU=refocus(uL(:,j),v,l(:,j),box_min,box_max,box_stp,lambda);
d2=size(trU,2); d2=ceil(d2/2);
figure, imshow(squeeze(abs(mean(trU(:,d2,:,:),4))),[])

parfor j=1:Nl
    j
    tic
    trU=refocus(uL(:,j),v,l(:,j),box_min,box_max,box_stp,lambda);
    rU(:,j)=trU(:);
    toc
    
    %d2=size(rU,2); d2=ceil(d2/2);
    %figure, imshow(squeeze(abs(rU(:,50,:,j))),[])
    %keyboard
end

dd=size(trU);
rU=reshape(rU,[dd,Nl]);
d2=size(rU,2); d2=ceil(d2/2);
figure, imshow(squeeze(abs(mean(rU(:,d2,:,:),4))),[])
