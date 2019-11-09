function [inliers, transf] = ransac(matches, blobs1, blobs2)
% This code is part of:
%
%   CMPSCI 670: Computer Vision, Fall 2019
%   University of Massachusetts, Amherst
%   Instructor: Subhransu Maji
%
%   Homework 4

% implement this

%% Get all the non-empty feature matches
[inputs, cols, bases] = find(matches);
% Reset max inliers counter
maxInliers = 0;
% Store the original indices of matches so that removing the 0 index entries doesn't change the actual indices
originalIndices = find(matches);
% Get the non-zero entries size
inputs_size = size(inputs,1);
if nargin< 4
    method = 'affine';
end

if strcmp(method,'affine')
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
    fprintf('Maximum Inliers %d and Msize is %d\n',maxInliers,inputs_size);
    transf = [inv(best_model(:,1:2)),-1*best_model(:,3)];
    inliers = originalIndices(inliersC);
    disp('Affine Transformation is:')
    disp(transf)
end

if strcmp(method, 'homography')
    % X co-ordinates of base points
    xBase = zeros(inputs_size,1);
    % Y co-ordinates of base points
    yBase = zeros(inputs_size,1);
    % X co-ordinates of transformed points
    xInput = zeros(inputs_size,1);
    % Y co-ordinates of transformed points
    yInput = zeros(inputs_size,1);
    
    % Blobs of Image 1 are Base Points
    for i= 1:inputs_size
        xBase(i,1) = blobs1(inputs(i),1);
        yBase(i,1) = blobs1(inputs(i),2);
    end
    
    % Blobs of Image 2 are Transformed Points
    for i= 1:inputs_size
        xInput(i,1) = blobs2(bases(i),1);
        yInput(i,1) = blobs2(bases(i),2);
    end
    
    
    for i = 1:10000
        
        randBase = zeros(4,2);
        randInput = zeros(4,2);
        
        % Picking up random base points for first step of RANSAC
        for j = 1:4
            index = ceil(rand*(inputs_size-1))+1;
            randBase(j,1) = xBase(index);
            randBase(j,2) = yBase(index);
            randInput(j,1) = xInput(index);
            randInput(j,2) = yInput(index);
            
        end
        
        h = computeHomography(randInput(:,1),randInput(:,2),randBase(:,1),randBase(:,2));
        
        [x,y] = applyHomography(inv(h),xBase,yBase);
        total = 0;
        for j = 1:numel(x)
            sigma = 100;
            error = sum(sum(([xInput(j) yInput(j)] - [x(j) y(j)]).^2));
            
            if  error < sigma
                total = total + 1;
            end
            
        end
        
        if total > maxInliers
            
            maxInliers = total;
            H = h;
            
        end
    end
    transf = H;
    inliers = [];
end
end

function [x2, y2] = applyHomography(H, x1,y1)
    p(1,:) = x1(:);
    p(2,:) = y1(:);
    p(3,:) = 1;

    x2 = H(1,:) * p;
    y2 = H(2,:) * p;
    thd  = H(3,:) * p;

    x2 = x2 ./ thd;
    y2 = y2 ./ thd;
end

function [H] = computeHomography(x1,y1,x2,y2)
    num = numel(x1);
    H = ones(3);
    eq = zeros(num*2, 8);
    sol = zeros(num*2,1);
    for i =1:num
        
        pts1 = [x1(i),y1(i),1];
        z = [0, 0, 0];
        rem = (-1 * [x2(i);y2(i)]) * [x1(i),y1(i)];
        sol((i-1)*2+1:(i-1)*2+2) = [x2(i); y2(i)];
        eq ((i-1)*2+1:(i-1)*2+2,:) = [[pts1, z;z, pts1] rem];
    end

H = eq\sol;

H(9) = 1;
H = reshape(H,3,3)';
end

function [dist]= calculateError(calPoints,FinputPoints)
    dist = sum((calPoints-FinputPoints).^2);
end

function [Transform] = computeAffine(basePoints,transPoints)
    baseP = zeros(6,6); % input points
    transP = zeros(6,1); % points after transformation of input points.
    for j = 1:3
        index = ceil(rand*(size(basePoints,1)-1))+1;
        % Prepare base matrix with points before transformation; put them in the position shown in affine transform matrix equation
        baseP(2*j-1,:) = [basePoints(index,1),basePoints(index,2),1,0,0,0];
        baseP(2*j,:) = [0,0,0,basePoints(index,1),basePoints(index,2),1];
        % Set transformed points matrice; put them in the position shown in affine transform matrix equation
        transP(2*j-1,:) = transPoints(index,1);
        transP(2*j,:) = transPoints(index,2);
    end
    % Matrix Division
    Transform = baseP\transP;
end
