% Calculate the right matrix
% points_of_left = dlmread('points_from_imL.txt');
% points_of_right = dlmread('points_from_imR.txt');

% Test the second pair of images
points_of_left = dlmread('points_from_imL_2.txt');
points_of_right = dlmread('points_from_imR_2.txt');

[num_of_points, num_of_cols] = size(points_of_left);

M = zeros(num_of_points, 9);

for row = 1:1:num_of_points
    left = [points_of_left(row, 1), points_of_left(row, 2), 1];
    right = [points_of_right(row, 1), points_of_right(row, 2), 1];
    temp = zeros(1, 9);
    count = 1;
    for i = 1:1:3
        for j = 1:1:3
            temp(count) = left(i) * right(j);
            count = count + 1;
        end
    end
    M(row, :) = temp; 
end

% disp('Big Matrix =')
% disp(M);

% Eight points method
% res = -1 * ones(8, 1);
% F = Eight_Points_Matrix \ res;

[U, S, V] = svd(M);
[min_val, min_index] = min(diag(S(:, :)));
f_m = V(:, min_index);

F = reshape(f_m, 3, 3);

disp('Fundamental Matrix =')
disp(F)

% pic_1 = imread('imL.png');
% Left image of second pair
pic_1 = imread('left_image.jpeg');
subplot(1, 2, 1), title('LEFT IMAGE'), imshow(pic_1), hold on;
% figure(1), imshow(pic_1), hold on;
for point = 1:1:num_of_points
    epipolar_line =  transpose(F) * transpose([points_of_right(point, :), 1.0]);
    a = epipolar_line(1);
    b = epipolar_line(2);
    c = epipolar_line(3);
    x = (0:0.001:384);
    y = - (a * x / b + c / b);
    scatter(points_of_left(point, 1), points_of_left(point, 2), 200, 'g.'), hold on;
    plot(x, y, 'LineWidth', 0.8), hold on;
end

% pic_2 = imread('imR.png');
% Right image of second pair
pic_2 = imread('right_image.jpeg');
subplot(1, 2, 2), title('RIGHT IMAGE'), imshow(pic_2), hold on;
for point = 1:1:num_of_points
    epipolar_line =  F * transpose([points_of_left(point, :), 1.0]);
    a = epipolar_line(1);
    b = epipolar_line(2);
    c = epipolar_line(3);
    x = (0:0.001:384);
    y = - (a * x / b + c / b);
    scatter(points_of_right(point, 1), points_of_right(point, 2), 200, 'g.'), hold on;
    plot(x, y, 'LineWidth', 0.8), hold on;
end

% Calculate the epipoles
[u, d] = eigs(transpose(F) * F);
uu = u(:, 1);
epipole = uu / uu(3);
subplot(1, 2, 1);
% Show the epipole
scatter(epipole(1), epipole(2), 500, 'b.');