function [pano] = MyPanorama()

%% YOUR CODE HERE.
% Must load images from ../Images/Input/
% Must return the finished panorama.

% actual program

images = [];
for ii = 1:9
    images = [images struct('im',double(imread(sprintf("../Images/TestSet2/%d.jpg", ii))))];
end


filt = fspecial('gaussian',5,0.5)/25;

H2 = stitchTwo(images(2).im, images(1).im, 300, 41, filt, 2, 1000, 6, 0.1, false);
H3 = stitchTwo(images(3).im, images(2).im, 300, 41, filt, 2, 1000, 6, 0.1, false);

warped2 = warpImage(images(2).im, H2);  
warped3 = warpImage(warpImage(images(3).im, H3), H2);  

[im1_height, im1_width, ~] = size(images(1).im);
[im2_height, im2_width, ~] = size(images(2).im);
[im3_height, im3_width, ~] = size(images(3).im);

b1 = blendTwo(warped2, images(1).im, H2, 1, im2_height, im2_width);
b2 = blendTwo(warped3, b1, [H3; H2], 2, im2_height, im2_width);

imshow(uint8(b2))

%}