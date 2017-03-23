imL = imread('imL.jpeg');

% Translate the color-scale image to grey-scale, prepare for the calculation of intensity
imL_grey = rgb2gray(imL);

imR = imread('imR.jpeg');

% Translate the color-scale image to grey-scale, prepare for the calculation of intensity
imR_grey = rgb2gray(imR);

[m, n] = size(imL_grey);

% Set up the matching block size
block_size = 51;
padding_size = (block_size - 1) / 2;

% (x1, y1) is the coordinate of selected point
points = dlmread('keyPoints.txt');

[row, col] = size(points);

subplot(1, 2, 1), title('LEFT IMAGE');
imshow(imL_grey), hold on;

subplot(1, 2, 2), title('RIGHT IMAGE');
imshow(imR_grey), hold on;

world_coordinates = zeros(row, 3);

for point = 1:1:row
    
    x = points(point, 1);
    y = points(point, 2);
    
    subplot(1, 2, 1);
    scatter(x, y, 'r.'), hold on;

    % This is a little trick, I assume this line (y = x + y1  0 < x < m) is the epipolar line
    epipolar_line = y;

    % Plot the rectangle of the point that need to be matched
    rectangle('Position', [x - padding_size  y - padding_size block_size block_size], 'EdgeColor','g'), hold on;

    % Select the matric that need to be matched
    selected_areas = zeros(block_size, block_size);
    for i = 1 : 1 : block_size
        for j = 1 : 1 : block_size
            selected_areas(i, j) = imL_grey(min(y + i - padding_size, n), min(x + j - padding_size, m));
        end
    end

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
                corresponding_matrix(i, j) = imR_grey(min(y + i - padding_size, m), min(d + j - padding_size, n));
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

    % Draw the two parallel lines to mark the sliding window
    subplot(1,2,2);
    x2 = 0 : 0.01 : n;
    y2 = x2 * 0 + epipolar_line;
    plot(x2, y2 - padding_size, 'r'), plot(x2, y2 + padding_size, 'r'), hold on;

    % Plot the matched point and its retangle
    plot(most_matched_point, y, 'r.'), hold on;
    rectangle('Position', [most_matched_point - padding_size  y - padding_size block_size block_size], 'EdgeColor','g'), hold on;
    % Calculate the disparity d * Z = f * B     f = 4.15 mm B = 20 mm
 
    % Through Assignment1, we know the pixel density, so that we know the relation between
    % pixel coordinate to mm coordinate
    pixel_unit_to_mm_unit = 1 / 833.3333;
    f = 4.15;
    B = 20;
    d = abs(most_matched_point - x) * pixel_unit_to_mm_unit;
%     d = abs(most_matched_point - x);
    Z = f * B / d;
    X = -(x * pixel_unit_to_mm_unit * Z) / f;
    Y = -(y * pixel_unit_to_mm_unit * Z) / f;
    world_coordinates(point, :) = [X, Y, Z];
end

% Store the world coordinates
dlmwrite('matched_points.txt', world_coordinates);

% Use the min value and max value of each coordinate to build the box that is similar to
% the object that I photographed before
% x_min = min(world_coordinates(:, 1));
% x_max = max(world_coordinates(:, 1));
% y_min = min(world_coordinates(:, 2));
% y_max = max(world_coordinates(:, 2));
% z_min = min(world_coordinates(:, 3));
% z_max = max(world_coordinates(:, 3));
% 
% figure(2);
% x=[x_min x_max x_max x_min x_min x_min x_min x_min x_min x_min x_min x_max x_max x_max x_max x_max x_max x_max x_min];
% y=[y_min y_min y_min y_min y_min y_min y_max y_max y_min y_max y_max y_max y_max y_min y_min y_max y_max y_max y_max];
% z=[z_min z_min z_max z_max z_min z_max z_max z_min z_min z_min z_max z_max z_min z_min z_max z_max z_min z_min z_min];
% plot3(x, y, z,'r'), hold on;

figure(2);
% 1 and 2
v = [world_coordinates(1, :); world_coordinates(2, :)];
plot3(v(:, 1),v(:, 2), v(:, 3), 'LineWidth', 1), hold on;
% 2 and 3
v = [world_coordinates(2, :); world_coordinates(3, :)];
plot3(v(:, 1),v(:, 2), v(:, 3), 'LineWidth', 1), hold on;
% 1 and 4
v = [world_coordinates(1, :); world_coordinates(4, :)];
plot3(v(:, 1),v(:, 2), v(:, 3), 'LineWidth', 1), hold on;
% 2 and 5
v = [world_coordinates(2, :); world_coordinates(5, :)];
plot3(v(:, 1),v(:, 2), v(:, 3), 'LineWidth', 1), hold on;
% 3 and 6
v = [world_coordinates(3, :); world_coordinates(6, :)];
plot3(v(:, 1),v(:, 2), v(:, 3), 'LineWidth', 1), hold on;
% 4 and 5
v = [world_coordinates(4, :); world_coordinates(5, :)];
plot3(v(:, 1),v(:, 2), v(:, 3), 'LineWidth', 1), hold on;
% 5 and 6
v = [world_coordinates(5, :); world_coordinates(6, :)];
plot3(v(:, 1),v(:, 2), v(:, 3), 'LineWidth', 1), hold on;