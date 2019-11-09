%% Compute the blob response
% Increase filter size, keep image the same

for i = 1:n
    sigma_next = sigma * k^(i-1); % Increament sigma for next level
    filt_size = 2*ceil(3*sigma_next)+1;  % Compute the filter kernel size in the same way as previously 
    LoG = sigma_next^2 * fspecial('log', filt_size, sigma_next);   % Making a custom LoG filter with trial and error
    filter_next = imfilter(im, LoG, 'same', 'replicate');  % Create a replica of the previous filter
    filter_next = filter_next .^ 2;     % Multiply it by itself
    scaleSpace(:,:,i) = filter_next;    % Store filter of each level     
end
