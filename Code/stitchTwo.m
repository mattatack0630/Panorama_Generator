% "fdbf" = "feature descriptor blurring filter"
function [homography] = stitchTwo(im1, im2, nbest, patchsize, fdbf, fmatchthresh, ransac_iter_limit, ransac_thresh, homog_fail_thresh, showOutput)
    %% find features (in each image)
    
    o_im1 = im1; o_im2 = im2;
    
    im1 = rgb2gray(uint8(im1)); im2 = rgb2gray(uint8(im2));
    
    C1 = cornermetric(im1);
    C2 = cornermetric(im2);
    
    %input("done finding cornermetrics; press enter");
    
    M1 = imregionalmax(C1); M2 = imregionalmax(C2);
    
    %input("got regionalmaxes; press enter");
    
    rawfeatures1 = getTruePositions(M1);
    rawfeatures2 = getTruePositions(M2);    
    
    % showing output for detected features
    if showOutput
        rawfeatposs = evertStructArray(rawfeatures1);
        figure
        imshow(uint8(o_im1));
        hold on 
        plot(rawfeatposs(:,1),rawfeatposs(:,2),'or');
        hold off
    end
    
    %input("got positions; press enter");
    
    %% extract N_best features (in each image)
    
    nbfeats1 = getNBest(C1, evertStructArray(rawfeatures1), nbest);
    nbfeats2 = getNBest(C2, evertStructArray(rawfeatures2), nbest);
    
    % showing output for features after ANMS
    if showOutput
        figure
        imshow(uint8(o_im1));
        hold on
        plot(nbfeats1(:,1),nbfeats1(:,2),'or');
        hold off
    end
    
    %input("got nbest features; press enter");
    
    %% create feature descriptors (for each image)
    
    padamt = (patchsize - 1)/2;
    padded1 = padarray(im1, [padamt padamt], "replicate");
    padded2 = padarray(im2, [padamt padamt], "replicate");
    
    %input("padded images; press enter");
    
    % the return values here should be structs, pairing feature position
    % structs with feature descriptors
    fds1 = makeFeatDescs(padded1, matrixToStructs(nbfeats1), patchsize, fdbf, 0);
    fds2 = makeFeatDescs(padded2, matrixToStructs(nbfeats2), patchsize, fdbf, 0);
    
    %input("created feature descriptors; press enter");
    
    %% find the best matches between features in each image
    
    featmatches = matchTheFeatures(fds1, fds2, fmatchthresh);
    
    if showOutput
        mf1 = [];
        mf2 = [];
        
        for mi = 1:length(featmatches)
            F1 = featmatches(mi).f1;
            F2 = featmatches(mi).f2;
            mf1 = [mf1 ; F1.x F1.y];
            mf2 = [mf2 ; F2.x F2.y];
        end
 
        showMatchedFeatures(im1, im2, mf1, mf2, 'montage');
    end
    
    %% estimate the best homography between the two images
    [homography,inliers] = findHomography(featmatches, ransac_iter_limit, ransac_thresh, homog_fail_thresh);
    
end

% return an array of structs representing the positions of each of the "1"
% values in a logical matrix
function [positions] = getTruePositions(logmat)
    [y,x] = size(logmat);
    positions = [];
    for xi = 1:x
        for yi = 1:y
            if logmat(yi,xi) == 1
                positions = [positions struct('x',xi,'y',yi)];
            end
        end
    end
end

% turn an array of structs into a struct of arrays
function [poss] = evertStructArray(positions)
    poss = [];
    num = length(positions);
    for pi = 1:num
        poss = [poss ; positions(pi).x positions(pi).y];
    end
end

function [sarray] = matrixToStructs(positions)
    sarray = [];
    [num ~] = size(positions);
    for pi = 1:num
        sarray = [sarray struct('x',positions(pi,1),'y',positions(pi,2))];
    end
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

% use ANMS to extract the N_best best corners. Takes a matrix of corner
% scores, an array of structs representing the positions of the detected
% corners, and a positive integer N_best. Returns an array of N_best
% positions, each of which is the position is a position of an optimal corner.

% REFACTOR:
% cscore doesn't change
% nbest doesn't change
% corners is now a n x 2 matrix representing corner positions
function [bestcorners] = getNBest(cscore, corners, nbest)
    
    [n ~] = size(corners);
    
    % this is unfortunately not a map but an array parallel to corners
    distToClosestSuperior = inf*ones(1,n);
    
    for ai = 1:n
        a = corners(ai,:);
        for bi = 1:n
            b = corners(bi,:);
            if cscore(b(2),b(1)) > cscore(a(2),a(1)) % checking if superior
                ED = (b(1) - a(1))^2 + (b(2) - a(2))^2; % computing distance
                
                if ED < distToClosestSuperior(ai)
                    distToClosestSuperior(ai) = ED;
                end
            end 
        end
    end
    
    %input("Done forming close-competitor map; press enter");
    
    bestcorners = [];
    % a quadratic approach - we can do much better but sorting
    % structs would require manually writing the sort
    
    while ((size(bestcorners) < nbest) == [1 1]) & n > 0
        % Find the max value of distToClosestSuperior and position of the
        % feature point parallel to it in the corners array. Then, we will
        % effectively ignore that maximum (while we search for the next
        % highest) by swapping it with the "last" element in the list and
        % decrementing the size of the list we're considering
        
        % 1: xpos
        % 2: ypos
        % 3: dist
        % 4: index        
        maxd = [corners(1,1) corners(1,2) distToClosestSuperior(1) 1];
        
        for ci = 2:n
            if maxd(3) < distToClosestSuperior(ci)
                maxd(3) = distToClosestSuperior(ci);
                maxd(1) = corners(ci,1);
                maxd(2) = corners(ci,2);
                maxd(4) = ci;
            end
        end
        
        % assign the new maximum distance for the subset of the array
        bestcorners = [bestcorners ; maxd(1) maxd(2)];
        
        % swap this current maximum and the last element, then decrement
        % the size we're considering
        corners(maxd(4),:) = corners(n,:);
        distToClosestSuperior(maxd(4)) = distToClosestSuperior(n);
        corners(n,:) = [maxd(1) maxd(2)];
        distToClosestSuperior(n) = maxd(3);
        n = n - 1;
        
    end
    
    %input("done extraction n_best from map; press enter");
end

% returns an array of structs pairing a feature position with its feature
% desciptor, based on the paddedim, a list of feature positions 'feats',
% a patch size, and a filter to apply to each patch before subsampling
function [fds] = makeFeatDescs(paddedim, feats, patchsize, filt, debug)
    fds = []; % initially empty array of feature descriptors
    [blay,blax] = size(paddedim);
    numfeats = length(feats); % get number of features
    fhps = (patchsize-1)/2;
    
    for fi = 1:numfeats
        cf = feats(fi); % go through each feature
        
        patch = zeros(patchsize,patchsize); % initialize empty patch
        
        for y = (cf.y - fhps):(cf.y + fhps) % iterate through pixels
            for x = (cf.x - fhps):(cf.x + fhps) % covered by the patch
    	        rely = y - cf.y + fhps + 1;
    	        relx = x - cf.x + fhps + 1;
                % and copy them to the patch matrix
                
                if debug == 1
                    [blay,blax]
                    y
                    x
                    rely
                    relx
                end
                
    	        patch(rely,relx) = paddedim(y + fhps, x + fhps);
            end
        end
    
        %blurredpatch = imfilter(patch,filt); % apply the filter ...
        %subsampled = imresize(blurredpatch,[8 8]); % subsample ...
        %fd = reshape(subsampled,[64 1]); % reshape ...
        %fd = zscore(fd); % and standardize.
        
        fd = zscore(reshape(imresize(imfilter(patch, filt), [8 8]), [64 1])); % (flattened code)
        
        fds = [fds struct('point',cf,'descriptor',fd)];
    end
end

% returns an array of feature matches, given two arrays of features with
% their descriptors, and a closeness threshold. fds1 and fds2 should be arrays
% of structs, having an attribute 'point' with a value that is a struct of
% a position (x,y attributes), and an attribute 'descriptor' with a value
% that is a 64-element column vector representing the feature descriptor.
function [matched_features] = matchTheFeatures(fds1, fds2, thresh)
    matched_features = []; % will be a struct array; the first element of 
    % the struct is the point in the first image; the second, the point in 
    % the second image
    numfeats1 = length(fds1);
    numfeats2 = length(fds2);
    
    % minimize the SSD for each match; using a threshold between the 1st
    % and 2nd best matches
    for fdi = 1:numfeats1
        fd1 = fds1(fdi);
        
        % initial minimum and "penminimum"
        best = struct('f1',struct('x',0,'y',0),'f2',struct('x',0,'y',0),'diff',inf);
        nextbest = struct('f1',struct('x',0,'y',0),'f2',struct('x',0,'y',0),'diff',inf);
        
        for fdj = 1:numfeats2
            fd2 = fds2(fdj);
            
            % find quality of match
            dista = fd2.descriptor - fd1.descriptor;
    	    dista = dot(dista, dista); % sum of square differences
    	    
            %disp(sprintf("attempting optimization"))
            % update maxima
            if dista < best.diff
                %disp(sprintf("found a next best: dista = %d <= best.diff = %d", dista, best.diff))
    	        nextbest = best;
    	        best.f1 = fd1.point; best.f2 = fd2.point; best.diff = dista;
            elseif dista < nextbest.diff
                %disp(sprintf("found a next best: dista = %d <= nextbest.diff = %d", dista, nextbest.diff))
    	        nextbest.f1 = fd1.point; nextbest.f2 = fd2.point; nextbest.diff = dista;
    	    end
        end
        
        % take match if it is sufficiently strong
        %disp(sprintf("final check:\n\tbest.diff = %d\n\tnextbest.diff = %d\n\t((nextbest.diff)/(best.diff)) = %d\n\tthresh = %d", ...
            %best.diff, nextbest.diff, ((best.diff)/(nextbest.diff)), thresh))
        if ((nextbest.diff)/(best.diff)) >= thresh
            matched_features = [matched_features struct('f1',best.f1,'f2',best.f2)];
        end
    end
end

% blend two images into one, now that im1 has been warped into the
% "coordinate space" of im2; the result of this warping is wim1. Our goal
% is to go from the coordinate space of the result of imwarp(im1) = wim1
% and the coordinate space of im2 to the coordinate space of the stitching
% field.
% H is the homography
% im1_height,im1_width are the respective dimensions of the unwarped image 1
function [blended] = blendTwo(wim1, im2, H, im1_height, im1_width)
    
    [wim1_height wim1_width ~] = size(wim1);
    [im2_height  im2_width  ~] = size(im2);
    
    % compute the wim1 corners in im2 coordinate space
    % does H need to be transposed first?
    [wim1_topleftx wim1_toplefty] = apply_homography(H,1,1);
    [wim1_bottomleftx wim1_bottomlefty] = apply_homography(H,1,im1_height);
    [wim1_toprightx wim1_toprighty] = apply_homography(H,im1_width,1);
    [wim1_bottomrightx wim1_bottomrighty] = apply_homography(H,im1_width,im1_height);
    
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

%{
function [idx] = saFind(structArray, query)
    idx = 1
    while idx <= length(structArray) & structArray(idx) ~= query
        idx = idx + 1;
    end
end
%}
function [H, best_inliners] = findHomography(feature_pairs, iter, thresh, failure_thresh)

    best_inliners = [];
    sz = length(feature_pairs);
    
    for i = 1:iter        
        [X, Y, x, y] = splitFeatures([ ...
            feature_pairs(randi(sz)), feature_pairs(randi(sz)), ...
            feature_pairs(randi(sz)), feature_pairs(randi(sz))]);
        
        homography = est_homography(X, Y, x, y); 
        inliners = [];
        
        for j = 1:sz
            p1 = feature_pairs(j).f1;
            p2 = feature_pairs(j).f2;
            [hx, hy] = apply_homography(homography, p1.x, p1.y);
            diff = [p2.x - hx, p2.y - hy];
            dist = diff(1)^2 + diff(2)^2;
            if dist < thresh
                inliners = [inliners; feature_pairs(j)];
            end    
        end
        
        if length(inliners) > length(best_inliners)
            best_inliners = inliners;
        end
        
        if length(best_inliners) >= 0.99 * sz
            break;
        end
    end
        
    disp(sprintf("RANSAC inlier/total ratio: %d",length(best_inliners)/sz));
    assert(length(best_inliners)/sz >= failure_thresh)
    
    [X, Y, x, y] = splitFeatures(best_inliners);
    H = est_homography(X, Y, x, y);    
end

function [X, Y, x, y] = splitFeatures(features)
    X = [];
    Y = [];
    x = [];
    y = [];
    
    for j = 1:length(features)
        fm = features(j);
        X = [X; fm.f2.x];
        Y = [Y; fm.f2.y];
        x = [x; fm.f1.x];
        y = [y; fm.f1.y];
    end
end