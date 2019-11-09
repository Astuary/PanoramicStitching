%% Get all the non-empty feature matches
[inputs, cols, bases] = find(matches);
% Reset max inliers counter
maxInliers = 0;
% Store the original indices of matches so that removing the 0 index entries doesn't change the actual indices
originalIndices = find(matches);
% Get the non-zero entries size
inputs_size = size(inputs,1);