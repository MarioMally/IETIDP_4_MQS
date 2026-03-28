function [nrbarr, condVec, dirHigh, dirLow, epsVec] = cnstrct_paper_example_config(dimCond, numLoops, type)
% Creates one layer of the inductor coil test example for the low frequency
% paper
%
%   nrbarr      array containing nurbs
%   condVec     vector containing information if the patch at the same index in
%               the nrbarr is a conductor or not
%
%   dimCond     dimension of the (square) conductor and also the air gaps
%   numLoops    number of conductor loops
if numLoops < 1
    error("not yet defined")
end

%% define values
thicknessCond = dimCond;
dimX = 3 + numLoops * 4;
dimY = 5 + numLoops * 4;

x0 = 0;
y0 = 0;

x1 = (dimX - 1)*thicknessCond;
y1 = ceil(dimY/2)*thicknessCond;

y2 = (dimY-1)*thicknessCond;

xend = dimX*thicknessCond;
yend = dimY*thicknessCond;

%% outer air layer
% builds the pattern of the outer air layer
clear outerLayer
outerLayer(1) = nrb4surf([x0 y0], [x0+thicknessCond x0+thicknessCond], [x0 y1], [x0+thicknessCond y1]);
outerLayer(2) = nrb4surf([x0 y1], [x0+thicknessCond y1], [x0 yend], [x0+thicknessCond y2]);
outerLayer(3) = nrb4surf([x0+thicknessCond y2], [x1 y2], [x0 yend], [xend yend]);
outerLayer(4) = nrb4surf([x1 y2], [xend yend], [x1 y1], [xend y1]);
outerLayer(5) = nrb4surf([x1 y1], [xend y1], [x1 y0+3*thicknessCond], [xend y0+3*thicknessCond]);
outerLayer(6) = nrb4surf([x1 y0+3*thicknessCond], [xend y0+3*thicknessCond], [x1 y0+2*thicknessCond], [xend y0+2*thicknessCond]);
outerLayer(7) = nrb4surf([x1 y0+2*thicknessCond], [xend y0+2*thicknessCond], [x1 y0+thicknessCond], [xend y0+thicknessCond]);
outerLayer(8) = nrb4surf([x1 y0+thicknessCond], [xend y0+thicknessCond], [x1 y0], [xend y0]);
outerLayer(9) = nrb4surf([x1 y0], [x0 y0], [x1 y0+thicknessCond], [x0+thicknessCond y0+thicknessCond]);
condVec = [0 0 0 0 0 0 1 0 0];

% outerLayer(1) = nrb4surf([x0 y0], [x0+a x0+a], [x0 y1], [x0+a y1]);
% outerLayer(2) = nrb4surf([x0 y1], [x0+a y1], [x0 yend], [x0+a y2]);
% outerLayer(3) = nrb4surf([x0+a y2], [x1 y2], [x0 yend], [xend yend]);
% outerLayer(4) = nrb4surf([x1 y2], [xend yend], [x1 y1], [xend y1]);
% outerLayer(5) = nrb4surf([x1 y1], [xend y1], [x1 y0+3*a], [xend y0+3*a]);


%% looped layers
% every iteration is conductor + inner layer of air
clear loopedLayersTmp
loopedLayers = [];

for i = 1:numLoops
    clear loopedLayersTmp

    % modify "constants"
    x0 = x0+thicknessCond;
    y0 = y0+thicknessCond;
    x1 = x1-thicknessCond;
    y2 = y2-thicknessCond;
    xend = xend-thicknessCond;
    yend = yend-thicknessCond;

    % simple conductor part
    loopedLayersTmp(1) = nrb4surf([x0 y0], [x0+thicknessCond y0+thicknessCond], [x0 y1], [x0+thicknessCond y1]);
    loopedLayersTmp(2) = nrb4surf([x0 y1], [x0+thicknessCond y1], [x0 yend], [x0+thicknessCond y2]);
    loopedLayersTmp(3) = nrb4surf([x0+thicknessCond y2], [x1 y2], [x0 yend], [xend yend]);
    loopedLayersTmp(4) = nrb4surf([x1 y2], [xend yend], [x1 y1], [xend y1]);
    loopedLayersTmp(5) = nrb4surf([x1 y1], [xend y1], [x1 y0+3*thicknessCond], [xend y0+2*thicknessCond]);

    % check for first iteration conductor part
    if i == 1
        loopedLayersTmp(6) = nrb4surf([xend y0], [x0 y0], [xend y0+thicknessCond], [x0+thicknessCond y0+thicknessCond]);
    else
        loopedLayersTmp(6) = nrb4surf([xend+2*thicknessCond y0], [x0 y0], [x1+2*thicknessCond y0+thicknessCond], [x0+thicknessCond y0+thicknessCond]);
    end

    % modify "constants" again for air part
    x0 = x0+thicknessCond;
    y0 = y0+thicknessCond;
    x1 = x1-thicknessCond;
    y2 = y2-thicknessCond;
    xend = xend-thicknessCond;
    yend = yend-thicknessCond;

    % simple air part
    loopedLayersTmp(7) = nrb4surf([x0 y0], [x0+thicknessCond y0+thicknessCond], [x0 y1], [x0+thicknessCond y1]);
    loopedLayersTmp(8) = nrb4surf([x0 y1], [x0+thicknessCond y1], [x0 yend], [x0+thicknessCond y2]);
    loopedLayersTmp(9) = nrb4surf([x0+thicknessCond y2], [x1 y2], [x0 yend], [xend yend]);
    loopedLayersTmp(10) = nrb4surf([x1 y2], [xend yend], [x1 y1], [xend y1]);

    % check for last iteration air part
    if i == numLoops
        loopedLayersTmp(11) = nrb4surf([x0+2*thicknessCond y1-thicknessCond], [x0+3*thicknessCond y1-thicknessCond], [x0+2*thicknessCond y1], [x0+3*thicknessCond y1]);
    else
        loopedLayersTmp(11) = nrb4surf([x1 y1], [xend y1], [x1 y0+3*thicknessCond], [xend y0+2*thicknessCond]);
    end

    % check for first iteration air part
    if i == 1
        loopedLayersTmp(12) = nrb4surf([xend+thicknessCond y0], [x0 y0], [xend+thicknessCond y0+thicknessCond], [x0+thicknessCond y0+thicknessCond]);
    else
        loopedLayersTmp(12) = nrb4surf([xend+2*thicknessCond y0], [x0 y0], [x1+2*thicknessCond y0+thicknessCond], [x0+thicknessCond y0+thicknessCond]);
    end

    % modify global arrays
    condVec = [condVec, 1 1 1 1 1 1 0 0 0 0 0 0];
    loopedLayers = [loopedLayers, loopedLayersTmp];
end

%% innermost layer

% modify "constants" one last time
x0 = x0+thicknessCond;
y0 = y0+thicknessCond;
x1 = x1-thicknessCond;
y2 = y2-thicknessCond;
xend = xend-thicknessCond;
yend = yend-thicknessCond;

clear innermostLayer
innermostLayer(1) = nrb4surf([xend+2*thicknessCond y0], [x0 y0], [x1+2*thicknessCond y0+thicknessCond], [x0+thicknessCond y0+thicknessCond]);
innermostLayer(2) = nrb4surf([x0 y0], [x0+thicknessCond y0+thicknessCond], [x0 y1], [x0+thicknessCond y1]);
innermostLayer(3) = nrb4surf([x0 y1], [x0+thicknessCond y1], [x0, y1+thicknessCond], [x0+thicknessCond,y1+thicknessCond]);

condVec = [condVec, 1 1 1];

%% combine
nrbsurfaces = [outerLayer, loopedLayers, innermostLayer];

%% extrude surfaces
nrbarr = [];
for i = 1:numel(nrbsurfaces)
    nrbarr = [nrbarr, nrbextrude(nrbsurfaces(i), [0, 0, thicknessCond])];
end
condVec_indu = condVec;

dirLow = 7;
dirHigh = nan;

%% return before adding air layers
if type == 0
    error("Boundary exitation is wrong for no air box outside!")
end

%% add air layer on bottom
condVec = [condVec, zeros(size(nrbsurfaces))];
condVec(end) = 1;
nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbarr(i), vectrans([0 0 -thicknessCond]))];
end
nrbarr = [nrbarr, nrbarr_tmp];

%% add air layer on top
condVec = [condVec, zeros(size(nrbsurfaces))];
nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbarr(i), vectrans([0 0 thicknessCond]))];
end
nrbarr = [nrbarr, nrbarr_tmp];

dirHigh = dirLow;
dirLow = numel(condVec) - 1/3*numel(condVec) ;
epsVec = ones(size(condVec));

%% return before adding second inductor
if type == 1
    return
end

if type == 2

%% add inductor layer
condVec = [condVec, condVec_indu];
nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbarr(i), vectrans([0 0 2*thicknessCond]))];
end
nrbarr = [nrbarr, nrbarr_tmp];

%% add air layer on top
condVec = [condVec, zeros(size(nrbsurfaces))];
condVec(end) = 1;
nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbarr(i), vectrans([0 0 3*thicknessCond]))];
end
nrbarr = [nrbarr, nrbarr_tmp];

dirHigh = 2*numel(nrbarr)/5;
end

%% return before adding capacitor
if type == 2
    return
end

fac = 25;

%% add first capacitor sheet
condVec_tmp = ones(size(nrbsurfaces));
condVec_tmp(1:9) = 0;
condVec = [condVec, condVec_tmp];
nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbtform(nrbarr(i), vecscale([1 1 1/fac])), vectrans([0 0 (-1-1/fac)*thicknessCond]))];
end
nrbarr = [nrbarr, nrbarr_tmp];


%% add air layer in between
epsVec = ones(size(condVec));
epsVec = [epsVec, 1000*ones(size(nrbsurfaces))];
condVec = [condVec, zeros(size(nrbsurfaces))];
nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbarr(i), vectrans([0 0 (-2-1/fac)*thicknessCond]))];
end
nrbarr = [nrbarr, nrbarr_tmp];

%% add second capacitor sheet
condVec_tmp = ones(size(nrbsurfaces));
condVec_tmp(1:9) = 0;
condVec = [condVec, condVec_tmp];
epsVec = [epsVec, ones(size(nrbsurfaces))];
nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbtform(nrbarr(i), vecscale([1 1 1/fac])), vectrans([0 0 (-2-2/fac)*thicknessCond]))];
end
nrbarr = [nrbarr, nrbarr_tmp];

%% add outermost air layer
condVec = [condVec, zeros(size(nrbsurfaces))];
condVec(end) = 1;
epsVec = [epsVec, ones(size(nrbsurfaces))];
nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbarr(i), vectrans([0 0 (-3-2/fac)*thicknessCond]))];
end
nrbarr = [nrbarr, nrbarr_tmp];

dirHigh = numel(condVec);

end
