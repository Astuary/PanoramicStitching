% Evaluation code for blob detection
% Your goal is to implement scale space blob detection using LoG 
%
% This code is part of:
%
%   CMPSCI 670: Computer Vision
%   University of Massachusetts, Amherst
%   Instructor: Subhransu Maji

% imageName = 'butterfly.jpg';
% imageName = 'einstein.jpg';
% imageName = 'fishes.jpg';
 imageName = 'sunflowers.jpg';

numBlobsToDraw = 1000;
imName = imageName(1:end-4);

dataDir = fullfile('..','data','blobs');
im = imread(fullfile(dataDir, imageName));


%% Implement the code to detect blobs here
blobs = detectBlobs(im); 

%% Draw blobs on the image
drawBlobs(im, blobs, numBlobsToDraw);
title('Blob detection');
