% Non-Maximum Suppression between Scales and Threshold
% For scales, compare the current one with the previous and next ones, choose the maximum of the three
for i = 1:n
    max_space(:,:,i) = max(max_space(:,:,max(i-1,1):min(i+1,n)),[],3);
end
% Zero Out All Positions that are not the Local Maxima of the score (if the Value if not Greater than all its Neighbors)
max_space = max_space .* (max_space == scaleSpace);