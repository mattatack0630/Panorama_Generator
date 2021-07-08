function [blended] = blendTwo(wim1, im2, H, hs, im1_height, im1_width)
    
    [wim1_height wim1_width ~] = size(wim1);
    [im2_height  im2_width  ~] = size(im2);
    
    % compute the wim1 corners in im2 coordinate space
    % does H need to be transposed first?
    wim1_topleftx = 0;
    wim1_toplefty = 0;
    wim1_bottomleftx = 0;
    wim1_bottomlefty = 0;
    wim1_toprightx = 0;
    wim1_toprighty = 0;
    wim1_bottomrightx = 0;
    wim1_bottomrighty = 0;

    for i = 1:hs
        h = H(i*3-2:(i)*3, :);
        [wim1_topleftx, wim1_toplefty] = apply_homography(h,1,1);
        [wim1_bottomleftx, wim1_bottomlefty] = apply_homography(h,1,im1_height);
        [wim1_toprightx, wim1_toprighty] = apply_homography(h,im1_width,1);
        [wim1_bottomrightx, wim1_bottomrighty] = apply_homography(h,im1_width,im1_height);
    end    
    
    % round each of the corner coordinates
    wim1_topleftx = round(wim1_topleftx);
    wim1_toplefty = round(wim1_toplefty);
    wim1_toprightx = round(wim1_toprightx);
    wim1_toprighty = round(wim1_toprighty);
    wim1_bottomleftx = round(wim1_bottomleftx);
    wim1_bottomlefty = round(wim1_bottomlefty);
    wim1_bottomrightx = round(wim1_bottomrightx);
    wim1_bottomrighty = round(wim1_bottomrighty);
    
    % the minimum x and y values in wim1, in im2's coordinate space
    wim1_xmin = min([wim1_topleftx wim1_bottomleftx]);
    wim1_ymin = min([wim1_toplefty wim1_toprighty]);
    
    % compute the stitching field extrema, in im2 coordinate space
    sf_xmin = min([wim1_xmin 1]);
    sf_xmax = max([wim1_toprightx wim1_bottomrightx im2_width]);
    sf_ymin = min([wim1_ymin 1]);
    sf_ymax = max([wim1_bottomlefty wim1_bottomrighty im2_height]);
    
    % initialize stitching field
    stitch_field = NaN * ones(sf_ymax - sf_ymin + 1, sf_xmax - sf_xmin + 1, 3);
        % will need to be converted to uint8 later
    
    % the transformation from the coordinate space of the warped image 1 to
    % the coordinates of the resulting blended image
    wim_to_sf = @(x,y) [x+wim1_xmin-sf_xmin y+wim1_ymin-sf_ymin];
        
    % copy wim1 to stitching field using the imwarp coords -> stitching space function
    for xi = 1:wim1_width
        for yi = 1:wim1_height
            dest_p = wim_to_sf(xi,yi); % compute destination for pixel in wim1
            stitch_field(dest_p(2),dest_p(1),:) = wim1(yi,xi,:); 
        end
    end
    
    % the transformation from the coordinate space of image 2 to the
    % coordinates of the resulting blended image
    im2_to_sf = @(x,y) [x+1-sf_xmin  y+1-sf_ymin];
    
    % copy im2 to stitching field using the simply offset transform
    % the transformation from the coordinate space of image 2 to the
    % coordinates of the resulting blended image
    
    % copy im2 to stitching field using the simply offset transform
    for xi = 1:im2_width
        for yi = 1:im2_height
            dest_p = im2_to_sf(xi,yi); % compute pixel destination
            if ~isnan(im2(yi,xi,:))
                if isnan(stitch_field(dest_p(2),dest_p(1),:)) % check if empty: then fill
                    stitch_field(dest_p(2),dest_p(1),:) = im2(yi,xi,:);
                else % otherwise average with the current occupant
                    %Overwriting test
                    stitch_field(dest_p(2),dest_p(1),:) = im2(yi,xi,:);
                    %stitch_field(dest_p(2),dest_p(1),:) = 0.5*(stitch_field(dest_p(2),dest_p(1),:) + im2(yi,xi,:));
                end
            end
        end
    end
    
    blended = stitch_field;
end