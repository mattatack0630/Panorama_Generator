function warpImg1 = warpImage(o_im1, homography)
    %% apply the homography and blend the images
    warpImg1 = o_im1;
    H = homography;
%    for H = homography
        tform = projective2d(transpose(H)); % format the homography correctly

        % we need to extract the NaN values (representing empty pixels) from
        % the original image 1 before warping, and then add these back in after
        % warping their positions.

        im1_empty_pixels = threedisnan(warpImg1); % find the empty pixels
        w_empty_pixels = imwarp(im1_empty_pixels, tform); % warp the empty pixels

        warpImg1 = double(uint8(warpImg1)); % cleanse the NaNs from the original
        warpImg1 = imwarp(warpImg1, tform, 'FillValues', NaN); % warp the cleansed image

        warpImg1 = emptyThePixels(warpImg1, w_empty_pixels); % incorporate the empty pixels back into the warped image
%    end
end


% this will set logmat(j,i) = true if im1(j,i,:) == [NaN NaN NaN]
%
% the purpose is to create a matrix which encodes the location of NaN in
% im1, which in turn encode "empty" pixels
function [logmat] = threedisnan(im1)
    [y x ~] = size(im1);
    logmat = zeros(y,x,'logical');
    for xi = 1:x
        for yi = 1:y
            if isnan(im1(yi,xi,:))
                logmat(yi,xi) = true;
            end
        end
    end
end

% this takes two matrices
% wim1 should be a 3-dimensional double matrix, and epls a 2-dimensional
% logical matrix with the same height and width as wim1
% this will make wim1(j,i,k) = NaN if epls(j,i) is true for all k
% epls = empty pixel locations
%
% the purpose of this function is the project is to allow use to preserve
% the locations of empty pixels in a source image, and warp the empty and
% nonempty pixels separately, and then add them together agains
function [warped] = emptyThePixels(wim1, epls)
    warped = wim1;
     
    [y x ~] = size(wim1);
    for xi = 1:x
        for yi = 1:y
            if epls(yi,xi) % if epls says there's an empty pixel here...
                warped(yi,xi,:) = [NaN NaN NaN]; % then put one there
            end
        end
    end
end
