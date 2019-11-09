% Try RANSAC for 'iterations' number of times
for itr = 1:iterations
    % Affine Transform for the current iteration
    T = computeAffine(inputPoints,basePoints);
    % Reshape the transform matrix in [3x3] for 
    % tranformed_points = tranform x base_points
    trans = [reshape(T,[3,2])';0,0,1];
    % Reshape base_points to [3x1] for
    % tranformed_points = tranform x base_points
    inputPt = [inputPoints';ones(1,inputs_size)];
    % Compute tranformed_points = tranform x base_points
    inputptF =  trans*inputPt;
    % Take the x-y co-ordinates of every tranformed point
    calPoints = inputptF(1:2,:);

    % Calculate Error between Expected Tranformed point co-ordinates and actual Transformed points
    dist = calculateError(calPoints,basePoints');
    % Only take the points whose distance is lower than threshold
    [rows,cols,error] = find(dist.*(dist < threshold));
    % Count the number of inliers
    inliersCount = size(cols,2);
    total_error = sum(error,2);

    % Find the best model by inlier count and error values
    if inliersCount > round(0.10*inputs_size) && inliersCount >= maxInliers && inliersCount ~= inputs_size

        best_model = reshape(T,[3,2])';
        if(inliersCount == maxInliers)
            maxIndex = find(inliersInfo(:,1) == maxInliers);
            if(total_error < inliersInfo(maxIndex,2))
                best_model = reshape(T,[3,2])';
            end
        end
        inliersC = cols';
        maxInliers= inliersCount;
        inliersInfo(itr,1)=maxInliers ;
        inliersInfo(itr,2)= total_error;
        transInfo(:,:,itr) = best_model;
    end
end
