function matches = computeMatches(f1,f2)
% COMPUTEMATCHES Match two sets of SIFT features f1 and f2
% This code is part of:
%
%   CMPSCI 670: Computer Vision, Fall 2019
%   University of Massachusetts, Amherst
%   Instructor: Subhransu Maji
%
%   Mini project 4

% Implement this
f1_size = size(f1, 1);
f2_size = size(f2, 1);

%% Computing Matches using SSD
matches = zeros(f1_size, 1);
% Start Stopwatch Timer
tic
if 0
    fprintf('Computing matching using SSD: \n');
    % For each feature from image 1
    for i = 1 : f1_size
       % Find its nearest match image 2 
       bestMatch = inf;
       for j = 1 : f2_size
            % SSD
            match = sum(sum((f1(i,:)-f2(j,:)).^2));
            % Store the lowest SSD feature
            if(match<bestMatch)
               matches(i) = j;
               bestMatch = match;
            end
       end
    end
% Stop Stopwatch Timer
toc
end

%% Computing Matches using ratio
if 1
    % For each Descriptor in the First Image, Select its Match to Second Image
    for i = 1 : f1_size
        % The most likely match
        bestMatch = inf;
        % The next most likely match after the best one
        secondMatch = inf;
        index = inf;
        for j = 1:f2_size
            % Compute SSD
            match = sum(sum((f1(i,:)-f2(j,:)).^2));
            % compare it with the best match
            if(match < bestMatch)
                secondMatch = bestMatch;
                bestMatch = match;
                index = j;
            % compare it with the second best match
            elseif(match < secondMatch && match ~= bestMatch)
                    secondMatch = match;
            end
        end
        % Now we have determined best and second best match
        
        % Computer Ratio
        ratio = bestMatch/secondMatch;
        % Compare it with a threshold 
        if(~isequal(index,inf) && ratio < 0.8)
            matches(i)= index;
        else
            matches(i) = 0;
        end
    end
end