%% Build far field covariance matrix

% Build the target area
boxTargetArea = boxArea( ...
    1 ,       ... wavelength
    10,       ... MFP (O.D. = 2)
    [-10,10], ... z
    [-10,10]  ... x
);

% views config, near field sensors, sensors located in a circle around the
% target area
viewsDirections = 0:1:359;
sensorsRadius = 20;
views = nearFieldSource(sensorsRadius, deg2rad(viewsDirections));

% lights config, lighting in some directions
lightsDirections = [0, 1, 4, 20];
lights = farFieldSource(deg2rad(lightsDirections),0);

% scatter config
% use Henyey-Greenstein scattering function
load('scatteringAmplitude.mat', 'theta', 'farField')
scatter = HGScatter(-0.5);

%% Sovle scmc with CBS
tic
CBSres = scmc(boxTargetArea,views,lights,scatter,1e3,'CBS',true,'parforIters',12);
toc

%% Sovle scmc without CBS
tic
NCBSres = scmc(boxTargetArea,views,lights,scatter,1e3,'CBS',false,'parforIters',12);
toc

%% Compare cov matrix
f = figure;
f.Position = [0,0,1200,700];
maxval = max(abs([CBSres.C(:);NCBSres.C(:)]));

subplot(2,4,1);
imagesc(viewsDirections,viewsDirections,abs(CBSres.C(:,:,1,1)),[0,maxval]);
xlabel('view[deg]');
ylabel({'CBS','view[deg]'});
title(['(',num2str(lightsDirections(1)),'\circ,',num2str(lightsDirections(1)),'\circ)']);

subplot(2,4,2);
imagesc(viewsDirections,viewsDirections,abs(CBSres.C(:,:,1,2)),[0,maxval]);
xlabel('view[deg]');
ylabel('view[deg]');
title(['(',num2str(lightsDirections(1)),'\circ,',num2str(lightsDirections(2)),'\circ)']);

subplot(2,4,3);
imagesc(viewsDirections,viewsDirections,abs(CBSres.C(:,:,1,3)),[0,maxval]);
xlabel('view[deg]');
ylabel('view[deg]');
title(['(',num2str(lightsDirections(1)),'\circ,',num2str(lightsDirections(3)),'\circ)']);

subplot(2,4,4);
imagesc(viewsDirections,viewsDirections,abs(CBSres.C(:,:,1,4)),[0,maxval]);
xlabel('view[deg]');
ylabel('view[deg]');
title(['(',num2str(lightsDirections(1)),'\circ,',num2str(lightsDirections(4)),'\circ)']);

subplot(2,4,5);
imagesc(viewsDirections,viewsDirections,abs(NCBSres.C(:,:,1,1)),[0,maxval]);
xlabel('view[deg]');
ylabel({'Without CBS','view[deg]'});

subplot(2,4,6);
imagesc(viewsDirections,viewsDirections,abs(NCBSres.C(:,:,1,2)),[0,maxval]);
xlabel('view[deg]');
ylabel('view[deg]');

subplot(2,4,7);
imagesc(viewsDirections,viewsDirections,abs(NCBSres.C(:,:,1,3)),[0,maxval]);
xlabel('view[deg]');
ylabel('view[deg]');

subplot(2,4,8);
imagesc(viewsDirections,viewsDirections,abs(NCBSres.C(:,:,1,4)),[0,maxval]);
xlabel('view[deg]');
ylabel('view[deg]');

%% Compare diagonals
figure

hold on
plot(viewsDirections,diag(abs(CBSres.C(:,:,1,1))), 'lineWidth', 2);
plot(viewsDirections,diag(abs(NCBSres.C(:,:,1,1))), 'lineWidth', 2);

xlabel('view[deg]');
ylabel('intensity');

legend('With CBS','Without CBS')