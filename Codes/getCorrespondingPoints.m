% first_image = imread('imL.png');
first_image = imread('left_image.jpeg');

imshow(first_image);

first_points = ginput();

% dlmwrite('points_from_imL.txt', first_points);
dlmwrite('points_from_imL_2.txt', first_points);

% second_image = imread('imR.png');
second_image = imread('right_image.jpeg');

imshow(second_image);

second_points = ginput();

% dlmwrite('points_from_imR.txt', second_points);
dlmwrite('points_from_imR_2.txt', second_points);