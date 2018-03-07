function [points, inliers, parameters, errors, dists] = fitting(fname,err_t, it, p_dist, dir)

fdir = 'D:\SZTAKI\Diploma\Pontfelho\dzsudzsabyte\';
file = strcat(fdir, fname);
% RANSAC parameters
noise_percent = 0.1;
iterations = it;% it; 200;
err_threshold = err_t;%0.005;%2e+19;%2.0e+19; %0.0005
err_dist = 10;
num_ind = 5;  %# points to fitting
alg = 0; % plane
mind = 1;
max_dist = 0.02;
nump = 1000;

points = importdata(file);

%r_idxs = datasample(1:length(points), nump);
%points = points(r_idxs(1:length(r_idxs)),:);

%points = prefilt_all(points, max_dist);
fout = 'plane_filt_frame_1600.xyz';
save_m(fout, points);

tStart = tic;
      
[inliers, outliers, parameters, errors, dists, sp, startp, pno_out] = RANSAC(points, iterations, err_threshold, num_ind, alg, mind, p_dist, 1);

tElapsed = toc(tStart)

fout = 'plane_in_frame_1600_kdtest.xyz';
save_m(fout, inliers);

%{
[inliers, outliers, parameters, errors, dists, sp, startp, pno_out] = RANSAC(outliers, iterations, err_threshold, num_ind, alg, mind, p_dist, 1);
fout = 'plane_in_frame_1600_it2_x0.xyz';
save_m(fout, inliers);

[inliers, outliers, parameters, errors, dists, sp, startp, pno_out] = RANSAC(outliers, iterations, err_threshold, num_ind, alg, mind, p_dist, 2);
fout = 'plane_in_frame_1600_it3_y0.xyz';
save_m(fout, inliers);

fout = 'plane_out_frame_1600_it3.xyz';
save_m(fout, outliers);
%}

%{
[p1 p2 p3 p4] = quaternion(points);
[inliers, parameters, errors, dists, sp] = RANSAC(p1, iterations, err_threshold, num_ind, alg, mind, p_dist);
[inliers2, parameters, errors, dists, sp] = RANSAC(p2, iterations, err_threshold, num_ind, alg, mind, p_dist);
[inliers3, parameters, errors, dists, sp] = RANSAC(p3, iterations, err_threshold, num_ind, alg, mind, p_dist);
[inliers4, parameters, errors, dists, sp] = RANSAC(p4, iterations, err_threshold, num_ind, alg, mind, p_dist);
%}
%[inliers, parameters, errors, dists, sp] = RANSAC([inliers; inliers2; inliers3; inliers4], iterations, err_threshold, num_ind, alg, mind, p_dist);


%{
fout = 'plane_in_frame_1600_p1.xyz';
save_m(fout, inliers);
fout = 'plane_in_frame_1600_p2.xyz';
save_m(fout, inliers2);
fout = 'plane_in_frame_1600_p3.xyz';
save_m(fout, inliers3);
fout = 'plane_in_frame_1600_p4.xyz';
save_m(fout, inliers4);
%}

%{
fout = 'plane_in_frame_3220_t.xyz';
save_m(fout, inliers);
fout = 'plane_frame_1600_selected_points.xyz';
save_m(fout, sp);
fout = 'plane_frame_1600_start_points.xyz';
save_m(fout, startp);
fout = 'plane_in_frame_start_points_orig.xyz';
save_m(fout, pno_out);
%}

function [inliers, outliers, parameters, errors, dists, sp, startp, pno_out] = RANSAC(points, iterations, err_threshold, num_ind, alg, mind, p_dist,dir)
    inliers=[];
    outliers = [];
    parameters = zeros(2, 3);
    indices = 1 : num_ind;
    errors = [];
    dists = [];
    sp = [];
    startp = [];
    pno_out = [];
    MdlKDT = KDTreeSearcher(points(:,1:3));
    while(iterations)
        %randomly/directly select + del
        [num_ind, points, p_act, pno] = select_points(points,p_dist,dir,MdlKDT);
        indices = 1 : num_ind;
        sp((length(sp) + 1), 1:3) = p_act;
        
        
        %if enough points to fit, find the parameters -> only 5 ponint -> unnesc.
        %{
        if (num_ind > 5) 
            num_ind = 5; 
            indices = 1 : num_ind;
        end
        %}
        
        if(num_ind == 5)
            switch alg
                case 0 %plane
                    [fitted, parameters_tmp] = fit_points_plane(points, indices);
                case 1 %sphere
                    [fitted, parameters_tmp] = fit_points_sphere(points, indices);
                case 2 %cilynder
                    [fitted, parameters_tmp] = fit_points_cilynder(points, indices);
                otherwise
                    inliers = -1;
                    outliers = -1;
                    parameters = -1;
                    errors = -1;
                    dists = -1;
                    return;
            end
        else
            fitted = false;
        end
        
        %
        if(fitted)          
            inliers_tmp = points(1:num_ind, :);
            %outliers_tmp = [];
            tmp_ind = num_ind;
            tmp_o_ind = 0;
            errors_tmp = [];
            dists_tmp = [];
            
            %count the inlier within a threshold in the remaining set of points
            %-> if fits increase the set -- error under the threshold
            for i = (num_ind+1) : length(points)
                err = error(points(i, 1:3), parameters_tmp);
                errors_tmp(length(errors_tmp) +1) = err;
            end
           
           %{
           emin = min(errors_tmp);
           emax = max(errors_tmp);
           lim = emin  + (emax-emin)*err_threshold;
           %}
           lim = err_threshold;
            
            for i = 1 : length(errors_tmp)
                if (errors_tmp(i) < lim) 
                    tmp_ind = tmp_ind + 1;
                    inliers_tmp(tmp_ind, :) = points(i + num_ind, :);
                else
                    %tmp_o_ind = tmp_o_ind + 1;
                    %outliers_tmp(tmp_o_ind,:) = points(i + num_ind, :);
                end
                
            end
            %if more inliers, then earlier -> store it
            
            %as more inlier az it possible
            %%{
            
            %meanp = mean(inliers_tmp);
            %[dists_tmp, d] = min_distance(meanp, inliers_tmp);
            
            if(length(inliers) < length(inliers_tmp)) % && d < mind)
                inliers = inliers_tmp;
                %outliers = outliers_tmp;
                parameters = parameters_tmp;
                errors = errors_tmp;
                pno_out = pno;
                if(length(inliers) < 5)
                    startp = inliers;
                else
                    startp = inliers(1:5, :);
                end
                %dists = dists_tmp;
            end
            %}
            
            %{
            idx = zeros(length(inliers_tmp));
            d = zeros(length(inliers_tmp));
            for i = 1:length(inliers_tmp)
                [idx(i),d(i)] = knnsearch(points(2:end,:),points(1,:), 'K', 1);
            end
            %}
            
            
            disp('fitted:');
            disp(iterations);
            disp('inliers:');
            disp(length(inliers));
        else
            disp('not fitted:');
            disp(iterations);
        end
        
        iterations = iterations -1;
    end
end

function [fitted, parameters] = fit_points_cylinder(points, indices)
    fitted = -1;
    parameters = -1;
end

function [fitted, parameters] = fit_points_sphere(points, indices)
    mean_p = mean(points);
    
end        

function [fitted, parameters] = fit_points_plane2(points, indices)
    
    mean_p = mean(points);
    covar_m = zeros(6,1);
    parameters = zeros(2, 3);
    fitted = false;
    num_p = length(points);
    
    for i = 1 : length(indices)
        idx = indices(i);
        diff = points(idx,:) - mean_p;
        covar_m(1) = (diff(1) * diff(1)) / num_p; %x*x
        covar_m(2) = (diff(1) * diff(2)) / num_p; %x*y
        covar_m(3) = (diff(1) * diff(3)) / num_p; %x*z
        covar_m(4) = (diff(2) * diff(2)) / num_p; %y*y
        covar_m(5) = (diff(2) * diff(3)) / num_p; %y*z
        covar_m(6) = (diff(3) * diff(3)) / num_p; %z*z
    end
    
    weighted_dir  = [0 0 0];
    axis_dir= [0 0 0];
    det_x = covar_m(4) * covar_m(6) - covar_m(5) * covar_m(5); % yy*zz - yz*yz
    det_y = covar_m(1) * covar_m(6) - covar_m(3) * covar_m(3); % xx*zz - xz*xz
    det_z = covar_m(1) * covar_m(4) - covar_m(4) * covar_m(4); % xx*yy - xy*xy
    
    %x dir
    axis_dir(1) = det_x;
    axis_dir(2) = (covar_m(3)*covar_m(5) - covar_m(2)*covar_m(6));
    axis_dir(3) = (covar_m(2)*covar_m(5) - covar_m(3)*covar_m(4));
    w = det_x * det_x;
    if(sum(weighted_dir .* axis_dir) < 0) 
        w = -w; 
    end
    weighted_dir = weighted_dir + axis_dir*w;
    
    %y dir
    axis_dir(1) = (covar_m(3)*covar_m(5) - covar_m(2)*covar_m(6));
    axis_dir(2) = det_y;
    axis_dir(3) = (covar_m(2)*covar_m(3) - covar_m(5)*covar_m(1));
    w = det_y * det_y;
    if(sum(weighted_dir .* axis_dir) < 0) 
        w = -w; 
    end
    weighted_dir = weighted_dir + axis_dir*w;
    
    %z dir
    axis_dir(1) = (covar_m(2)*covar_m(5) - covar_m(4)*covar_m(3));
    axis_dir(2) = (covar_m(2)*covar_m(3) - covar_m(1)*covar_m(5));
    axis_dir(3) = det_z;
    w = det_z * det_z;
    if(sum(weighted_dir .* axis_dir) < 0) 
        w = -w; 
    end
    weighted_dir = weighted_dir + axis_dir*w;
    
    fitted = true;
    parameters(1,:) = mean_p(1:3) ;
    parameters(2,:) = weighted_dir ;
    
end 

function [fitted, parameters] = fit_points_plane(points, indices)
    points_in(:,:) = points(indices(:), 1:3);
    mean_p = mean(points_in); %points
    covar_m = zeros(6,1);
    parameters = zeros(2, 3);
    fitted = false;
    
    for i = 1 : length(indices)
        idx = indices(i);
        diff = points(idx,:) - mean_p;
        covar_m(1) = covar_m(1) + diff(1) * diff(1); %x*x
        covar_m(2) = covar_m(2) + diff(1) * diff(2); %x*y
        covar_m(3) = covar_m(3) + diff(1) * diff(3); %x*z
        covar_m(4) = covar_m(4) + diff(2) * diff(2); %y*y
        covar_m(5) = covar_m(5) + diff(2) * diff(3); %y*z
        covar_m(6) = covar_m(6) + diff(3) * diff(3); %z*z
    end
    
    
    det_x = covar_m(4) * covar_m(6) - covar_m(5) * covar_m(5); % yy*zz - yz*yz
    det_y = covar_m(1) * covar_m(6) - covar_m(3) * covar_m(3); % xx*zz - xz*xz
    det_z = covar_m(1) * covar_m(4) - covar_m(4) * covar_m(4); % xx*yy - xy*xy
    %det = max([det_x det_y, det_z]);
    det = max([det_x det_y, det_z]);
    
    if(det > 0)
        %inv_det = 1/det;
        parameters(1,:) = mean_p(1:3) ;
        %parameters(2,1) = (covar_m(4)*covar_m(3) - covar_m(2)*covar_m(5))*inv_det;
        %parameters(2,2) = (covar_m(1)*covar_m(5) - covar_m(2)*covar_m(3))*inv_det;
        %parameters(2,3) = -1;
        if(det == det_x)
            parameters(2,1) = det_x;
            parameters(2,2) = (covar_m(3)*covar_m(5) - covar_m(2)*covar_m(6));
            parameters(2,3) = (covar_m(2)*covar_m(5) - covar_m(3)*covar_m(4));
        elseif(det == det_y)
            parameters(2,1) = (covar_m(3)*covar_m(5) - covar_m(2)*covar_m(6));
            parameters(2,2) = det_y;
            parameters(2,3) = (covar_m(2)*covar_m(3) - covar_m(5)*covar_m(1));
        else
            parameters(2,1) = (covar_m(2)*covar_m(5) - covar_m(4)*covar_m(3));
            parameters(2,2) = (covar_m(2)*covar_m(3) - covar_m(1)*covar_m(5));
           parameters(2,3) = det_z;      
        end  
        
        parameters(2,:) = parameters(2,:) / norm(parameters(2,:));
        
        fitted = true;
    end
end 

function err = error(point, parameters)
    err = (sum((point - parameters(1,:)) .* parameters(2,:) ))^2;
end

function save_m(fname, my_matrix)
    fid = fopen(fname,'wt');
    for i = 1:size(my_matrix,1)
        fprintf(fid,'%g ',my_matrix(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end

function [d, md] = min_distance(point, set)
    d = zeros(length(set), 1);
    for i = 1:length(set)
        d(i) = sqrt( sumsqr (set(i,:) - point) );
    end
    %md = mean(mink(d, 200));
    md=mink(d, 200);
end


function array = mink(array, num)
    %%{
    array =  sort(array);
    if(length(array) > num)
            array = array(num);
            %array = array(1:num);
    else
        array = array(length(array));
    end
    %}
    
end

function [num_ind, points, p_act, pnearsorig] = select_points(points,p_dist, dir, MdlKDT)
    if (dir == 4)
        [num_ind, points, p_act, pnearsorig] = select_points_rand(points);
    else
        [num_ind, points, p_act, pnearsorig] = select_points_dir(points,p_dist, dir, MdlKDT) ;  
    end
end


function [num_ind, points, p_act, pnearsorig] = select_points_rand(points)
    
    %%{
    rows = randperm(size(points,1));
    points = points(rows,:);
    num_ind = 4;
    p_act = points(1,:);
    pnearsorig = [];
    %}
end   

function [num_ind, points, p_act, pnearsorig] = select_points_dir(points, p_dist, dir, MdlKDT)
    %%{
    % Select a random point
    idx = datasample(1:length(points), 1);
    point = points(idx,:);
    p_act = point;
    
    % Find the 5 suggested point in the near
    %{
    p_nears_x0 = get_p_nears_x0(point,p_dist);
    p_nears_y0 = get_p_nears_y0(point,p_dist);
    p_nears = p_nears_x0;
    pnearsorig = p_nears_x0;
    %}
    switch dir
        case 1
            p_nears = get_p_nears_x0(point,p_dist);
        case 2
            p_nears = get_p_nears_y0(point,p_dist);
        case 3
            p_nears = get_p_nears_y0(point,p_dist);
    end
    
     pnearsorig = p_nears;
    
    % Find the closest real point to the 5 point
    idxs = zeros(5,1);
    idxs(1) = idx;

    for i = 1:4
       %[idxs(i+1),d] =  find_closep(points(:,1:3),p_nears(i,1:3));
       
       
       
       

       %ez volt a jó
       %[idxs(i+1),d] = knnsearch(points(:,1:3),p_nears(i,1:3));
      
       [idxs(i+1), d] = knnsearch(MdlKDT,p_nears(i,1:3));
       
      
       
       if (d > 1)
           idxs(i+1) = idxs(1);
       end
    end
    
    % Delete duplicatios
    idxs = unique(idxs);
    num_ind = length(idxs);
    
    % Reorder the points
    % - Reorder the indicies
    idxs_out = 1:length(points);
    for i=1:length(idxs)
        idxs_out(idxs(i) - (i-1)) = [];
    end
    idxs = [idxs' idxs_out];
    % - Reorder the points
    points_new = points;
    for i=1:length(idxs)
    points_new(i,:) = points(idxs(i),:);
    end
    points = points_new;
    
    %}
    
    
    %rows = randperm(size(points,1));
    %points = points(rows,:);
    %while()
    %end
    
    %{
    idx1 = datasample(1:length(points), 1);
    points([1, idx1]) = points([idx1, 1]);
    idx2 = datasample(2:length(points), 1);
    points([2, idx2]) = points([idx2, 2]);
    idx3 = datasample(3:length(points), 1);
    points([3, idx3]) = points([idx3, 3]);
    idx4 = datasample(4:length(points), 1);
    points([4, idx4]) = points([idx4, 4]);
    
    while( abs(points(1,3)-points(2,3)) < 0.1 )
        idx2 = datasample(3:length(points), 1);
        points([2, idx2]) = points([idx2, 2]);
    end
    
    while( abs(points(1,3)-points(3,3)) < 0.1 )% || abs(points(2,3)-points(3,3)) < 0.1 )
        idx3 = datasample(4:length(points), 1);
        points([3, idx3]) = points([idx3, 3]);
    end
    
     while( abs(points(1,3)-points(4,3)) < 0.1)%  || abs(points(2,3)-points(4,3)) < 0.1 || abs(points(3,3)-points(4,3)) < 0.1 )
        idx4 = datasample(5:length(points), 1);
        points([4, idx4]) = points([idx4, 4]);
     end
    %}

    %{
    [idx,d] = knnsearch(points(2:end,:),points(1,:), 'K', num_ind);
    idxs_out = 1:length(points)';
    for i=1:length(idx)
    idxs_out(idx(i)) = [];
    end
    idxs = [idx idxs_out];
    points_new = points;
    for i=1:length(idxs)
    points_new(i) = points(idxs(i));
    end
    points = points_new;
    %}

    %points_tmp = points;
    %points_tmp((1:num_ind), :) = points(d(:),:);
    
    
    %{
    rows = randperm(size(points,1));
    points = points(rows,:);
    num_ind = 4;
    p_act = points(1,:);
    pnearsorig = [];
    %}
end


function p_nears = get_p_nears_x0(point,p_dist)  
   
    %4 point version
    %%{
    p_diffs = [0 -p_dist -p_dist ;
               0 -p_dist p_dist ; 
               0 p_dist -p_dist ;
               0 p_dist p_dist  ];
    point_4 = repmat(point, 4, 1);
    p_nears = point_4 + p_diffs;
    %}
end

function p_nears = get_p_nears_y0(point,p_dist)
   
    %4 point version
    %%{
    p_diffs = [-p_dist 0 -p_dist;
                -p_dist 0 p_dist ; 
                p_dist 0 -p_dist ;
                p_dist 0 p_dist];
    point_4 = repmat(point, 4, 1);
    p_nears = point_4 + p_diffs;
    %}
end

function p_nears = get_p_nears_z0(point,p_dist)
   
    %4 point version
    %%{
    p_diffs = [-p_dist -p_dist 0;
                -p_dist p_dist 0; 
                p_dist -p_dist 0;
                p_dist p_dist 0];
    point_4 = repmat(point, 4, 1);
    p_nears = point_4 + p_diffs;
    %}
end


function p_nears = get_p_nears_8(point,p_dist)
    
    %8 point version
    %%{
    p_diffs = [-p_dist -p_dist -p_dist; -p_dist -p_dist p_dist; 
                -p_dist p_dist -p_dist; -p_dist p_dist p_dist; 
                p_dist -p_dist -p_dist; p_dist -p_dist p_dist; 
                p_dist p_dist -p_dist; p_dist p_dist p_dist ];
    point_8 = repmat(point, 8, 1);
    p_nears = point_8 + p_diffs;
    %}

end


function [p1 p2 p3 p4] = quaternion(points)
    p1 = [];
    p2 = [];
    p3 = [];
    p4 = [];
    
    for i = 1 : length(points);
        if(points(i,1) > 0)
            if(points(i, 2) > 0)
                p1(length(p1) + 1, 1:3 ) = points(i,:);
            else
                p2(length(p2) + 1, 1:3 ) = points(i,:);
            end
        else
            if(points(i, 2) > 0)
                p3(length(p3) + 1, 1:3 ) = points(i,:);
            else
                p4(length(p4) + 1, 1:3 ) = points(i,:);
            end
        end
    end
end


function new_points = prefilt_all(points, max_dist)
    new_points = [];
    for i = 1 : length(points)
        new_points(i, 1:3) = prefilt(points(i,:), points, max_dist);
        disp(i);
    end
end

function point = prefilt(point, points, max_dist)
    near_p = [];
    for i = 1 : length(points)
        if(pdist([point;points(i,:)]) < max_dist)
            near_p(length(near_p) + 1, 1:3) = points(i,:);
        end
    end
    
    if(length(near_p) > 1 )
        point = mean(near_p);
    end
    
end

function [idx_close,d] = find_closep(points, point)
    dist =  [];
    for i = 1 : length(points)
        dist(i) = sqrt( (points(i,1)-point(1))^2 +  (points(i,2)-point(2))^2 + (points(i,3)-point(3))^2 );
    end
    d = min(dist);
    idx_close = find(dist == d); 
end

end