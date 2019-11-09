function [x2, y2] = applyHomography(H, x1,y1)
    p(1,:) = x1(:);
    p(2,:) = y1(:);
    p(3,:) = 1; % Homogenous co-ordinate

    % Compute Transformed Point
    x2 = H(1,:) * p;
    y2 = H(2,:) * p;
    thd  = H(3,:) * p;

    x2 = x2 ./ thd;
    y2 = y2 ./ thd;
end