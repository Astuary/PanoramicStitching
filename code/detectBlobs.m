function blobs = detectBlobs(im, param)
% DETECTBLOBS detects blobs in an image
%   BLOBS = DETECTBLOBSCALEIMAGE(IM, PARAM) detects multi-scale blobs in IM.
%   The method uses the Laplacian of Gaussian filter to find blobs across
%   scale space.
% 
% Input:
%   IM - input image
%   PARAM - struct containing the following fields
%       PARAM.SIGMA - sigma of the LoG filter (smallest blob desired)
%       PARAM.INTERVAL - number of intervals in an octave
%       PARAM.THRESHOLD - threshold for blob detection
%       PARAM.DISPLAY - if true then then shows intermediate results
%
% Ouput:
%   BLOBS - n x 4 array with blob in each row in (x, y, radius, score)
%

% Convert image to grayscale and convert it to double [0 1].
if size(im, 3) > 1
    im = rgb2gray(im);
end
if ~isfloat(im)
    im = im2double(im);
end

% If param is not specified use default
if nargin < 2
    param.sigma = 2;
    param.interval = 12;
    param.threshold = 1e-2;
    param.display = false;
end

% dummy blob
% blobs = [size(im, 2)/2, size(im, 1)/2, 100, 1.0];


%% Implement these:

%% Compute the scale space representation
sigma       = 1.6;      % Initial Scale 
k           = sqrt(2);  % Factor by which the scale is mulitplied each time  
sigma_final = power(k,16);      % Last Scale
% Dynamically decide the interation levels from the first and last scale values and multiplication factor
n           = ceil((log(sigma_final) - log(sigma))/log(k)); % Iterations of Laplacian scale space
[h, w]      = size(im); 
scaleSpace  = zeros(h, w, n);   % [h,w] - Dimensions of Image, n - Number of Levels in Scale Space

% Generate the Laplacian of Gaussian for the First Scale Level
filt_size   = 2 * ceil(3*sigma) + 1; 
LoG         = fspecial('log', filt_size, sigma);

%% Compute the blob response
% Increase filter size, keep image the same
if 0
    tic
    for i = 1:n
        sigma_next = sigma * k^(i-1); % Increament sigma for next level
        filt_size = 2*ceil(3*sigma_next)+1;  % Compute the filter kernel size in the same way as previously 
        LoG       =  sigma_next^2 * fspecial('log', filt_size, sigma_next);   % Making a custom LoG filter with trial and error
        filter_next = imfilter(im, LoG, 'same', 'replicate');  % Create a replica of the previous filter
        filter_next = filter_next .^ 2;     % Multiply it by itself
        scaleSpace(:,:,i) = filter_next;    % Store filter of each level     
    end
    toc
end

% Faster version: keep the filter size, downsample the image
if 1 
    tic
    imRes = im; % Create a copy of the image
    for i = 1:n
        imFiltered = imfilter(imRes, LoG, 'same', 'replicate'); % Apply the Laplacian of Gaussian filter here, the filter is already created in the above section
        imFiltered = imFiltered .^ 2; % Square of the LoG  
        scaleSpace(:,:,i) = imresize(imFiltered, size(im), 'bicubic'); % Resize the image to original size with bicubic interpolation
        if i < n        
            imRes = imresize(im, 1/(k^i), 'bicubic'); % Decrease the image size for next iteration
        end
    end
    toc
end

%% Perform NMS (spatial and scale space)
% Perform Non-Maximum Suppression for each Scale-Space Slice
suppression_size = 3; % To make it Maximum Filter
max_space = zeros(h, w, n); % NMS output for each level
for i = 1:n
    max_space(:,:,i) = ordfilt2(scaleSpace(:,:,i), suppression_size^2, ones(suppression_size));
end

% Non-Maximum Suppression between Scales and Threshold
% For scales, compare the current one with the previous and next ones, choose the maximum of the three
for i = 1:n
    max_space(:,:,i) = max(max_space(:,:,max(i-1,1):min(i+1,n)),[],3);
end
% Zero Out All Positions that are not the Local Maxima of the score (if the Value if not Greater than all its Neighbors)
max_space = max_space .* (max_space == scaleSpace);

x = [];   % X - axis of blob's center 
y = [];   % Y - axis of blob's center
radius = []; % Radius of the blob
score = [];     

% For each level
for i=1:n
    % Set a Threshold on the Squared Laplacian Response above which to report Region Detections
    [rows, cols, value] = find(max_space(:,:,i).*(max_space(:,:,i) >= 0.01));
    numBlobs = length(rows);
    % Create Lists of each tuple fields separately 
    if(numBlobs > 0)
        radii =  sigma * k^(i-1) * sqrt(2); 
        radii = repmat(radii, numBlobs, 1);
        x = [x; rows];
        y = [y; cols];
        score = [score;value];
        radius = [radius; radii];
    end
end

blobs = [y, x, radius, score];
