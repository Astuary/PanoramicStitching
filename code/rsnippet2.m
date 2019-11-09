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