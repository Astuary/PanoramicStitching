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