x = [];   % X - axis of blob's center 
y = [];   % Y - axis of blob's center
radius = []; % Radius of the blob
score = [];     

% For each level
for i=1:n
    % Set a Threshold on the Squared Laplacian Response above which to report Region Detections
    [rows, cols, value] = find(max_space(:,:,i).*(max_space(:,:,i) >= 0.01));
    numBlobs = length(rows);
    % Create Lists of each tuple fields separately 
    if(numBlobs > 0)
        radii =  sigma * k^(i-1) * sqrt(2); 
        radii = repmat(radii, numBlobs, 1);
        x = [x; rows];
        y = [y; cols];
        score = [score;value];
        radius = [radius; radii];
    end
end

blobs = [y, x, radius, score];