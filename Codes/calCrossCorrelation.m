imL = imread('imL.jpeg');

% Translate the color-scale image to grey-scale, prepare for the calculation of intensity
imL_gray = rgb2gray(imL);

imR = imread('imR.jpeg');
% Resize the image
% imR = imresize(imR, size);
% Translate the color-scale image to grey-scale, prepare for the calculation of intensity
imR_gray = rgb2gray(imR);

figure(1);
subplot(1,2,1);
imshow(imL_gray), hold on;
% [m, n] = size(imL_grey);

% (x1, y1) is the coordinate of selected point
[x1, y1] = ginput(1);
scatter(x1, y1, 'r.'), hold on;

% This is a little trick, I assume this line (y = x + y1  0 < x < m) is the epipolar line
epipolar_line = y1;

% Round the float value to int value
x1 = int32(x1);
y1 = int32(y1);

% Set up the matching block size
block_size = 21;
padding_size = (block_size - 1) / 2;

% Plot the rectangle of the point that need to be matched
rectangle('Position', [x1 - padding_size  y1 - padding_size block_size block_size], 'EdgeColor','g'), hold on;

% Select the matric that need to be matched
selected_areas = zeros(block_size, block_size);
for i = 1 : 1 : block_size
    for j = 1 : 1 : block_size
        selected_areas(i, j) = imL_gray(y1 + i - padding_size, x1 + j - padding_size);
    end
end

% Draw the two parallel lines to mark the sliding window
subplot(1,2,2);
imshow(imR_gray), hold on;
[m, n] = size(imR_gray);
x2 = 0 : 0.01 : n;
y2 = x2 * 0 + epipolar_line;
plot(x2, y2 - padding_size, 'r'), plot(x2, y2 + padding_size, 'r'), hold on;

% calculate the cross correlation and find the most matched point
selected_areas_zeromean = (selected_areas - mean(selected_areas(:))) ./ var(selected_areas(:));
normalized_crosscorrelation = zeros(n - block_size + 1);
ncc = -1;
max_ncc = -1;
most_matched_point = 0;

% Calcalate the cross correlation along the epipolar line
for d = padding_size : 1 : n - padding_size - 1
    corresponding_matrix = zeros(block_size, block_size);
    for i = 1 : 1 : block_size
        for j = 1 : 1 : block_size
            corresponding_matrix(i, j) = imR_gray(y1 + i - padding_size, d + j - padding_size);
        end
    end
    corresponding_matrix_zeromean = (corresponding_matrix - mean(corresponding_matrix(:))) ./ var(selected_areas(:));
    normalized = sqrt(sum(dot(selected_areas_zeromean, selected_areas_zeromean)) * sum(dot(corresponding_matrix_zeromean, corresponding_matrix_zeromean))); 
    ncc = sum(dot(selected_areas_zeromean, corresponding_matrix_zeromean)) / normalized;
    normalized_crosscorrelation(d - padding_size + 1) = ncc;
    if ncc > max_ncc
        max_ncc = ncc;
        most_matched_point = d;
    end
end

% Plot the matched point and its retangle
rectangle('Position', [most_matched_point - padding_size  y1 - padding_size block_size block_size], 'EdgeColor','g'), hold on;
plot(most_matched_point, y1, 'r.'), hold on;

% Plot the curve of the cross correlation
figure(2), title('Cross Correlation');
d = (padding_size : 1 : n - padding_size - 1);
plot(d, normalized_crosscorrelation, 'r');
disp(max_ncc);
disp(most_matched_point);

% Plot the stereo map and disparity map
figure(3);
subplot(1,2,1);
imshow(stereoAnaglyph(imL_gray,imR_gray));
subplot(1,2,2);
disparityRange = [-6 10];
disparityMap = disparity(imL_gray, imR_gray, 'BlockSize', 15, 'DisparityRange', disparityRange);
imshow(disparityMap);

% Calculate the disparity d * Z = f * B     f = 4.15 mm B = 20 mm
f = 4.15;
B = 20;
d = abs(most_matched_point - x1);
Z = f * B / d;
disp(Z);
X = - (x1 * Z) / f;
Y = - (y1 * Z) / f;
world_coordinate = [X, Y, Z];
disp(world_coordinate);