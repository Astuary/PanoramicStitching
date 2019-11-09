%% Computing Matches using ratio

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