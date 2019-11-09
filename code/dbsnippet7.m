imRes = im; % Create a copy of the image
for i = 1:n
    imFiltered = imfilter(imRes, LoG, 'same', 'replicate'); % Apply the Laplacian of Gaussian filter here, the filter is already created in the above section
    imFiltered = imFiltered .^ 2; % Square of the LoG  
    scaleSpace(:,:,i) = imresize(imFiltered, size(im), 'bicubic'); % Resize the image to original size with bicubic interpolation
    if i < n        
        imRes = imresize(im, 1/(k^i), 'bicubic'); % Decrease the image size for next iteration
    end
end