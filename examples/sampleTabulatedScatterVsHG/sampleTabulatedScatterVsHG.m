%% Sample Tabulated Scattering vs HG
% In this example, we will show how to build tabulated scatter in 2D and 3D
% and compare the sampling from the tabulated scattering scattering in
% multiple scattering case

% first build the target area for 2D and 3D
boxTargetArea2D = boxArea( ...
    1 ,         ... wavelength
    200,        ... MFP
    [-300,300], ... z
    [-300,300]  ... x
);

boxTargetArea3D = boxArea( ...
    1 ,         ... wavelength
    200,        ... MFP
    [-300,300], ... z
    [-300,300], ... x
    [-300,300]  ... y
);

% and the lighting and views
viewDirections = 0:1:360; % in deg
views = farFieldSource(deg2rad(viewDirections),0);
lights = farFieldSource(0,0); % light in 0 deg direction

% the g parameter we comapre with
gParam = 0.7;

%% Build 2D tabulated HG function

% the 2D direction vector MUST being with 0 and end with 2*pi
directions2D = (0:1e-4:1) * 2 * pi; 
hg2Damplitude = sqrt(evaluateHG(directions2D, gParam, 0, 2));

% solve for both tabulated and HG
tic
hgRes = scmc(boxTargetArea2D, views, lights, HGScatter(gParam), 1e4, ...
    'parforIters', 12);
toc

tic
tabRes = scmc(boxTargetArea2D, views, lights, ...
    tabulatedAmplitudeScatter(directions2D,hg2Damplitude), 1e4, ...
    'parforIters', 12);
toc

% plot the intensity of both results
figure
polarplot(deg2rad(viewDirections),diag(abs(hgRes.C)));
hold on
polarplot(deg2rad(viewDirections),diag(abs(tabRes.C)));

legend('HG','Tabulated HG');
title('2D HG plot')

%% Build 3D tabulated HG function

% the 3D direction vector MUST being with 0 and end with pi, theta is the
% elevation direction
cosThetaVals3D = (0:1e-4:1) * pi; 
hg3Damplitude = sqrt(evaluateHG(cosThetaVals3D, gParam, 0, 3));

% solve for both tabulated and HG
tic
hgRes = scmc(boxTargetArea3D, views, lights, HGScatter(gParam), 1e4, ...
    'parforIters', 12);
toc

tic
tabRes = scmc(boxTargetArea3D, views, lights, ...
    tabulatedAmplitudeScatter(cosThetaVals3D,hg3Damplitude), 1e4, ...
    'parforIters', 12);
toc

% plot the intensity of both results
figure
polarplot(deg2rad(viewDirections),diag(abs(hgRes.C)));
hold on
polarplot(deg2rad(viewDirections),diag(abs(tabRes.C)));

legend('HG','Tabulated HG');
title('3D HG plot')