% Convert image to grayscale and convert it to double [0 1].
if size(im, 3) > 1
    im = rgb2gray(im);
end
if ~isfloat(im)
    im = im2double(im);
end
