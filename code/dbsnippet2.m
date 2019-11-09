%% Compute the scale space representation
sigma = 1.6;  % Initial Scale 
k = sqrt(2);  % Factor by which the scale is mulitplied each time  
sigma_final = power(k,16); % Last Scale
% Dynamically decide the interation levels from the first and last scale values and multiplication factor
n = ceil((log(sigma_final) - log(sigma))/log(k)); % Iterations of Laplacian scale space
[h, w] = size(im); 
scaleSpace = zeros(h, w, n); % [h,w] - Dimensions of Image, n - Number of Levels in Scale Space

% Generate the Laplacian of Gaussian for the First Scale Level
filt_size = 2 * ceil(3*sigma) + 1; 
LoG = fspecial('log', filt_size, sigma);
