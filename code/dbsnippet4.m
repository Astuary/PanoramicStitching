% Perform Non-Maximum Suppression for each Scale-Space Slice
suppression_size = 3; % To make it Maximum Filter
max_space = zeros(h, w, n); % NMS output for each level
for i = 1:n
    max_space(:,:,i) = ordfilt2(scaleSpace(:,:,i), suppression_size^2, ones(suppression_size));
end