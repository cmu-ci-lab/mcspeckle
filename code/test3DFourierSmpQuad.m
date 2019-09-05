
lambda=1;
box_min=[-50; -50;  0];
box_max=[ 50;  50;  100];
box_min=[-4000; -4000;  0];
box_max=[ 4000;  4000;  100];
box_min=[-1000;-1000;0]; box_max=[1000;1000;30];%3000



sigt=1/2000;
sigt=1/25;
%theta_l=deg2rad([-35:5:35]);
theta_l=deg2rad([-20:20:20]);

theta_max=deg2rad(40);
theta_stp=deg2rad(0.2);
theta_max=0.5; 
theta_max=0.02; 
theta_stp=theta_max/200;
theta_l=[-theta_max/2:theta_stp*50:theta_max/2];


%IMPORTANT: this is a uniform grid of sin(theta) not of theta
%This is important to have a uniform grid for the Fourier transform
theta_v=[-theta_max:theta_stp:theta_max];

Nl0=length(theta_l);
Nv0=length(theta_v);
Nv=Nv0^2;
Nl=Nl0^2;



[vx,vy]=ndgrid(theta_v,theta_v);
v=-[vx(:)';vy(:)'; sqrt(1-vx(:).^2-vy(:).^2)'];  %Note again: we assume the grid is a uniform grid of sin(theta) so no sin applyied  


[lx,ly]=ndgrid(theta_l,theta_l);
l=[lx(:)';ly(:)'; sqrt(1-lx(:).^2-ly(:).^2)'];  


v(3,:)=-v(3,:); l(3,:)=-l(3,:);



doCBS=1; smpFlg=1;
maxItr=10^4*4;
maxItr=10^3
%maxItr=3;
%maxItr=100;
%maxItr=50;
%maxItr=2;
j=1

rng(520)
for j=1:1
tic

[uLf]=MCfieldFourierQuad( sigt, 1, box_min,box_max, l, theta_max,theta_stp,maxItr,lambda,doCBS,smpFlg,[],median(sign(v(3,:))));

toc
end
uLf=sum(uLf,4);

rng(520)

box_min=[-1007;-1007;0];
box_max=[ 1007; 1007;30];


for j=1:1
tic
uL(:,:,j)=MCfield( sigt, 1, box_min,box_max, l, v,1,1,maxItr,lambda,doCBS,smpFlg,1,0);
toc
end
uL=sum(uL,3);
figure, imshow(real(uLf(:,:,1)),[])
figure, imshow(real(reshape(uL(:,1),Nv0,Nv0)),[])
su=reshape(sum(uL,2),Nv0,Nv0);
suf=sum(uLf,3);
csuf=fftCorr(suf);
csu=fftCorr(su);



figure, imshow(abs(csu),[])
figure, imshow(abs(csuf),[])


csuf=fftCorr(uLf(:,:,13),uLf(:,:,15));
csu=fftCorr(reshape(uL(:,13),Nv0,Nv0),reshape(uL(:,15),Nv0,Nv0));
figure, imshow(abs(csu),[])
figure, imshow(abs(csuf),[])
csuf=fftCorr(uLf(:,:,13),uLf(:,:,1));
csu=fftCorr(reshape(uL(:,13),Nv0,Nv0),reshape(uL(:,1),Nv0,Nv0));
figure, imshow(abs(csu),[])
figure, imshow(abs(csuf),[])


return
