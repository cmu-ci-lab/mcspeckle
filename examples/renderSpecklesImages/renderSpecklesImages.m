%% Render speckles images
% Show how to use the rendering algorithm, and demonstrate the memory
% effect

% Build the target area
boxTargetArea = boxArea( ...
    1 ,             ... wavelength
    200,            ... MFP
    [0,2000],       ... z
    [-50000,50000], ... x
    [-50000,50000]  ... y
);

% g = 0.9 for strong ME. g = 0 will give weak ME.
scatter = HGScatter(0.9, 0.8);

% Light directions number, each direction is a different image
lightsNum = 5;

% Delta beween each direction (in degrees)
lightDelta = 0.0025;

lights = (0:1:(lightsNum-1)) * lightDelta;

% Create 101x101 image
pixNum = 101;

% For each lighting direction, we will set up 20 sensors
% The more sensorsToLight is bigger, the bigger speckles are, but less
% speckles will appear.
% It is also the numbet of pixels each speckle will move between lightings.
sensorsToLight = 20;

viewsVector = (-((pixNum-1)/2):1:((pixNum-1)/2)) * ...
    (lightDelta/sensorsToLight);

% Make views image
[X,Y] = meshgrid(viewsVector);

% The directions will be auto-normalized in the code
views = [sind(X(:).');sind(Y(:).');ones(1,numel(X))];

%% Run the code and get the rendered image
% Be careful not render the correlation matrix, since it has huge memory
% consuption

tic;
mulres = scmc(boxTargetArea, ...
    farFieldSource(views), ...
    farFieldSource(deg2rad(lights),0), ...
    scatter,1e3, ...
    'renderCov', false, ... Avoid from making the correlation matrix
    'CBS',false, ... Since we forward, the CBS is not necessary
    'uniformFirstScatter', false, ... For more efficient sampling
    'renderField', true, ... In order to actually render the image
    'parforIters', 12);
toc;

% put back the resolt into an image
u = reshape(mulres.field,pixNum,pixNum,lightsNum);

%% Plot the images
f = figure;
f.Position = [0,0,1400,250];

maxVal = max(abs(u(:)));

for lNum = 1:1:lightsNum
    subplot(1,lightsNum,lNum)
    imagesc(abs(u(:,:,lNum)),[0,maxVal]);
    
    colormap hot
    xticks([]);
    yticks([]);
    
    title([num2str(lights(lNum)),'\circ'])
end