figure(1);
image_left = imread('imL.jpeg');
image_left_grey = rgb2gray(image_left);
imshow(image_left_grey);
keyPoints = ginput();
dlmwrite('keyPoints.txt', keyPoints);
