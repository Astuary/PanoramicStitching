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