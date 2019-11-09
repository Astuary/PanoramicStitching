% This structure stores current point co-ordinates 
inputPoints = zeros(inputs_size,2);
% This structure stores the points whose co-ordinates are getting transformed
basePoints = zeros(inputs_size,2);
% Number of points whose distance from the line is less than some constant threshold
threshold = 4;
iterations = 35;
% Store the best affine transformation (with the least error)
best_model = zeros(2,3);
% Store the points who are inlier
inliersInfo = zeros(inputs_size,2);
% Store the affine transformation for a particular iteration
transInfo = zeros(2,3,inputs_size);

% Store all features and locations in inputPoints for 1st image and in basePoints for 2nd image
for ind= 1:inputs_size
    basePoints(ind,1) = blobs1(bases(ind),1);
    basePoints(ind,2) = blobs1(bases(ind),2);
    inputPoints(ind,1) = blobs2(inputs(ind),1);
    inputPoints(ind,2) = blobs2(inputs(ind),2);
end