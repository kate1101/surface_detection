function fitting_planes(file_name_in)

% Parameters in
iterations_fix = 100;
error_threshold = 0.01;

pointset_min_resolution = 0; %!!
pointset_min_size = 0; %!!
synthetic_point_distance = 0.5;


% Parameters out
fdir_in = 'D:\SZTAKI\Diploma\Pontfelho\dzsudzsabyte\';
fdir_out = 'D:\SZTAKI\Diploma\Code\Tests\T7\';
file_in = fullfile(fdir_in,file_name_in); 
file_name = strsplit(file_name_in,'.');
file_name_out = ''; %strcat(fdir_out, file_name(1), '_plane_1.', file_name(2));
idx_plane = 0; %1;
file_parameters = strcat(fdir_out,  file_name(1), '_parameters.txt');
all_parametrers= [];

% Load pointcloud
points = importdata(file_in);
points = points(:,1:3);
max_point_num = length(points)/2;
KDTS = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fitting_spheres;
%find_planes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nested functions

function find_planes
    % Finds and saves the planes of a pointcloud, and saves the outliers
    pointset_resolution = 1; %!!
    pointset_size = 1; %!!
    %while(pointset_resolution > pointset_min_resolution && pointset_size > pointset_min_size)
    for j = 1:6
        disp('planes:');
        disp(j);
        plane = RANSAC_pane;
        %separate_planes(plane);
        get_size_and_resolution(plane);
        set_new_plane_name;
        save_points(file_name_out, plane);
    end
    
    set_outlier_name;
    save_points(file_name_out, points);
    
    function separate_planes(plane)
        % Find more surface in the same plane, if its possible, and save them
        % into file
        get_size_and_resolution(plane);
        %+num_planes!!
    end
    
    function get_size_and_resolution(plane)
        pointset_resolution = 01; %!!
        pointset_size = 0; %!!
    end
end

function plane = RANSAC_pane
    % Detects the plane with the most inliers
    %points_to_search = [];
    plane = [];
    idxs_in = [];
    parameters = [];
    idxs_in_best = [];   
    idxs_out_best = [];
    error = 0;
    iterations = iterations_fix;
    % If there is too mutch input point, select max_point_num random
    % points, this makes the nearest neighbour point searh faster
    if( length(points) > max_point_num)
        row_idxs = datasample(1:length(points), max_point_num);
        %points_to_search = points(row_idxs(1:length(row_idxs)),:);
    else
        row_idxs = 1:length(points);
        %points_to_search = points;
    end
    
    KDTS = KDTreeSearcher(points(row_idxs(1:length(row_idxs)),:));
    
    % While iterations is true...
    while(iterations)
        disp('iterations:');
        disp(iterations);
        
        %... select points...
        if (idx_plane == 5)
            idxs_in = select_ground_points;
        else
            idxs_in = select_points;
        end
        %... fit a plane...
        [fitted, parameters] = plane_fitting( points(idxs_in(:),:) );
        %... if fitted, calc the errors and save the in/outliers
        if(fitted)
            idxs_in = []; % 1:length(points);
            for i = 1 : length(points)%1 : length(points)%length(points) : -1 : 1
                calc_error(points(i,:));
                if (error < error_threshold)
                   idxs_in(length(idxs_in) + 1) = i;
                end
            end
            %... if this is the best result, save 
            if( length(idxs_in) > length(idxs_in_best) )
                idxs_in_best = idxs_in;
                disp('inliers:');
                disp(length(idxs_in_best));
            end
        end
        iterations = iterations -1;
    end
    
    % Set the new plane
    plane = points(idxs_in_best(:),:);
    % Delete the plane points from tghe pointcloud
    
    if (idxs_in_best(1) ~= 1)
        idxs_out_best = [idxs_out_best 1 : (idxs_in_best(1)-1)];
    end
    
    for i = 1 : (length(idxs_in_best)-1)
        idxs_out_best = [idxs_out_best (idxs_in_best(i)+1) : (idxs_in_best(i+1)-1) ];
    end
    
    if (idxs_in_best(length(idxs_in_best)) ~= length(points))
        idxs_out_best = [idxs_out_best (idxs_in_best(length(idxs_in_best))+1) : length(points) ];
    end
    
    points = points(idxs_out_best(:),:);
    
    function idxs_in = select_ground_points
        idxs_in = []; 
        idxidx = datasample(1:length(row_idxs), 1);
        point = points(row_idxs(idxidx), :);
        p_z0 = get_p_nears_z0(point, synthetic_point_distance);
        [idxs_in, d_z0] = knnsearch(KDTS,p_z0);       
    end
    
    function idxs_in = select_points
        idxs_in = []; 
        idxidx = datasample(1:length(row_idxs), 1);
        point = points(row_idxs(idxidx), :);
        [p_x0, p_y0, p_z0] = get_synthetic_points(point, synthetic_point_distance);
        [idxs_near_x0, d_x0] = knnsearch(KDTS,p_x0);
        [idxs_near_y0, d_y0] = knnsearch(KDTS,p_y0);
        
        if (min(d_x0) < 0.05)
            d_x0 = 0;
        end
        
        if (min(d_y0) < 0.05)
            d_y0 = 0;
        end       
        
        if(sum(d_x0) > sum(d_y0))
            idxs_in = idxs_near_x0;
        else
            idxs_in = idxs_near_y0;
        end        
    end
       
    %whole version  
    %{
    function idxs_in = select_points_3dir
        idxs_in = []; 
        idxidx = datasample(1:length(row_idxs), 1);
        point = points(row_idxs(idxidx), :);
        [p_x0, p_y0, p_z0] = get_synthetic_points(point, synthetic_point_distance);
        [idxs_near_x0, d_x0] = knnsearch(KDTS,p_x0);
        [idxs_near_y0, d_y0] = knnsearch(KDTS,p_y0);
        [idxs_near_z0, d_z0] = knnsearch(KDTS,p_z0);
        
        if (min(d_x0) < 0.05)
            d_x0 = 0;
        end
        
        if (min(d_y0) < 0.05)
            d_y0 = 0;
        end
        
        if (min(d_z0) < 0.05)
            d_z0 = 0;
        end
        
        d_x0 = sum(d_x0);
        d_y0 = sum(d_y0);
        d_z0 = sum(d_z0);
        d = max([d_x0 d_y0 d_z0]);
        
        if(d == d_x0)
            idxs_in = idxs_near_x0;
        elseif(d == d_y0)
            idxs_in = idxs_near_y0;
        else
            idxs_in = idxs_near_z0;
        end
        
    end
    %}
    function calc_error(point)
        error = (sum((point - parameters(1,:)) .* parameters(2,:) ))^2;
    end
end

function [p_x0, p_y0, p_z0] = get_synthetic_points(point, p_dist)
    p_x0 = get_p_nears_x0(point,p_dist);  
    p_y0 = get_p_nears_y0(point,p_dist);  
    p_z0 = get_p_nears_z0(point,p_dist);  
end

function p_nears_x0 = get_p_nears_x0(point,p_dist)  
   
    %4 point version
    %%{
    p_diffs = [0 -p_dist -p_dist ;
               0 -p_dist p_dist ; 
               0 p_dist -p_dist ;
               0 p_dist p_dist  ];
    point_4 = repmat(point, 4, 1);
    p_nears_x0 = point_4 + p_diffs;
    %}
end

function p_nears_y0 = get_p_nears_y0(point,p_dist)
   
    %4 point version
    %%{
    p_diffs = [-p_dist 0 -p_dist;
                -p_dist 0 p_dist ; 
                p_dist 0 -p_dist ;
                p_dist 0 p_dist];
    point_4 = repmat(point, 4, 1);
    p_nears_y0 = point_4 + p_diffs;
    %}
end

function p_nears_z0 = get_p_nears_z0(point,p_dist)
   
    %4 point version
    %%{
    p_diffs = [-p_dist -p_dist 0;
                -p_dist p_dist 0; 
                p_dist -p_dist 0;
                p_dist p_dist 0];
    point_4 = repmat(point, 4, 1);
    p_nears_z0 = point_4 + p_diffs;
    %}
end

function [fitted, parameters] = plane_fitting(selected_points)
    mean_p = mean(selected_points); %points
    covar_m = zeros(6,1);
    parameters = zeros(2, 3);
    fitted = false;
    
    for i = 1 : length(selected_points)
        diff = selected_points(i,:) - mean_p;
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
    det = max([det_x det_y, det_z]);
    
    if(det > 0)
        parameters(1,:) = mean_p;
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
    
function save_points(file_name_out, actual_points)
    % Save a matrix into file
    fid = fopen(char(file_name_out),'wt');
    for i = 1:size(actual_points,1)
        fprintf(fid,'%g ',actual_points(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end

function set_new_plane_name
    idx_plane = idx_plane + 1;
    file_name_out = strcat(fdir_out,  file_name(1), '_plane_' , num2str(idx_plane), '.', file_name(2));
end

function set_outlier_name
    file_name_out = strcat(fdir_out, file_name(1), '_out.', file_name(2));
end

function save_parameters(parameters, file_parameters)
    fid = fopen(char(file_parameters),'wt');
    for i = 1:size(actual_points,1)
        fprintf(''
        fprintf(fid,'%g ',actual_points(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end

end