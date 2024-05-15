clear all;
close all;
clc
FNT_SZ = 20;
img = imread('../data/PalazzoTe.jpg');
figure(1);imshow(img);title('Original image');
%% convert rgb to grayscale and histogram equalization
img_gray=rgb2gray(img);
img_gray=imrotate(img_gray, -90);
figure(2); imshow(img_gray);
img_gray_eq=histeq(img_gray);
figure(3); imshow(img_gray_eq);
figure(4);subplot(2,1,1); histogram(img_gray); title('Grayscale histogram')
subplot(2,1,2)
histogram(img_gray_eq);title('Grayscale equalized image histogram')
%% canny edge detection
edgs = edge(img_gray_eq, 'canny',[.1 .2],10);
% the first vector is the bounding of the hysteresis thresholding and the 
% number is sigma (standard deviation of the gaussian used for the smooothing)
figure(5); imshow(edgs);title('Canny edge detection')
%% Harris Corner detection
corners = detectHarrisFeatures(img_gray_eq,'MinQuality', 0.13);
figure(6);imshow(img_gray_eq); hold on;
plot(corners);hold off











