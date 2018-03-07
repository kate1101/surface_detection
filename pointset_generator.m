function [s, point_cloud_p, c, new_s, new_c, inliers, parameters, errors] = pointset_generator

%Output params
s = [];
point_cloud_p = [];
c = [];
new_s = []; 
new_c = [];
inliers = [];
parameters = [];
errors = [];

%# of surface points (inliers)
num_in_1 = 500;
num_in_2 = 4000;
%# of outlier points 
num_out = 500;

% Parameters of the sphere 
t = [10, 23, 38];
R = 100;

% Parameters of the plane
P = [1, 2, 1];
A = [3, 5, 1];
B = [2, 3, 3];
u_max = 100;
v_max = 60;

% Parameters of the cylinder
t2 = [-7, 10, 1];
R2 = 75;
u2_max = 200;

% RANSAC parameters
noise_percent = 0.1;
iterations = 1000; 
err_threshold = 0.005; %2e+19;%2.0e+19; %0.0005
num_ind = 4; %# points to fitting
alg = 0; % plane

s1 = get_sphere(t, R, num_in_2);
s2 = get_noisy( s1 , noise_percent );
s3 = mutilation(s1);
s4 = mutilation(s2);

save_m('sphere_10_23_38_100.xyz', s1);
save_m('sphere_10_23_38_100_n.xyz', s2);
save_m('sphere_10_23_38_100_m.xyz', s3);
save_m('sphere_10_23_38_100_n_m.xyz', s4);

%{
s = get_noisy( get_sphere(t, R, num_in_2) , noise_percent );
p = get_noisy( get_plane(P, A, B, u_max, v_max, num_in_1) , noise_percent );
c = get_noisy( get_cylinder(R2, t2, u2_max, num_in_2) , noise_percent );
new_s = mutilation(s);
new_c = mutilation(c);

save_m('sphere.xyz', s);
save_m('plane.xyz', p);
save_m('cylinder.xyz', c);
save_m('new_sphere.xyz', new_s);
save_m('new_cylinder.xyz', new_c);

point_cloud_p = outliers(p, num_out);
point_cloud_c = outliers(new_c, num_out);
point_cloud_s = outliers(new_s, num_out);

save_m('plane_noisy.xyz', point_cloud_p);
save_m('cylinder_noisy.xyz', point_cloud_c);
save_m('sphere_noisy.xyz', point_cloud_s);

[inliers, parameters, errors] = RANSAC(point_cloud_p, iterations, err_threshold, num_ind, alg);
save_m('plane_in.xyz', inliers);
%}

%------------------------Nested functions----------------------------------

% Generate plane
function p = get_plane(p, a, b, u_max, v_max, num)
% [x] = [p] + [a]u + [b]v
% u, v : {R}

u = u_max * rand(num, 1);
v = v_max * rand(num, 1);
x = p(1) + (a(1) .* u) + (b(1) .* v); 
y = p(2) + (a(2) .* u) + (b(2) .* v); 
z = p(3) + (a(3) .* u) + (b(3) .* v); 

rgb = repmat(255, num, 3);
p = [x , y, z, rgb];

end

% Generate sphere
function s = get_sphere(t, r, num)
% u: [0;2pi]
% v: [0;pi]

u = 2 * pi * rand(num, 1);
v = pi * rand(num, 1);
x = r .* sin(v) .* cos(u) + t(1); 
y = r .* sin(v) .* sin(u) + t(2);
z = r .* cos(v) + t(3);

rgb = repmat(255, num, 3);
s = [x , y, z, rgb];
end

function c = get_cylinder(r, t, u_max, num)
% u : R
% v : [0,2pi]

u = u_max * rand(num, 1);
v = 2 * pi * rand(num, 1);
x = r .* cos(v) + t(1); 
y = r .* sin(v) + t(2);
z = u + t(3);

rgb = repmat(255, num, 3);
c = [x , y, z, rgb];
end

function new_surface = get_noisy(surface, noise_percent)
    
    new_surface = surface;
    max_val = max(surface(1:3));
    min_val = min(surface(1:3));
    diff = max_val - min_val;
    noise = noise_percent * max(diff);
    
    for i = 1 : length(surface)
       new_surface(i, 1) = surface(i, 1) + (rand - 0.5) * noise;
       new_surface(i, 2) = surface(i, 2) + (rand - 0.5) * noise;
       new_surface(i, 3) = surface(i, 3) + (rand - 0.5) * noise;
    end
end

function new_surface = mutilation(surface)
    new_surface = 0;
    max_val = max(surface);
    min_val = min(surface);
    margins = min_val + (max_val - min_val)/2;
    j = 1;
    for i = 1: length(surface)
        if sum(margins <= surface(i,:)) == 6
            new_surface(j,1:6) = surface(i,:);
            j = j+1;
        end
    end
    
end

function save_m(fname, my_matrix)
    fid = fopen(fname,'wt');
    for i = 1:size(my_matrix,1)
        fprintf(fid,'%g ',my_matrix(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end

function point_cloud = outliers(surface, out_num)
    max_val = max(max(surface(:,1:3)));
    min_val = min(min(surface(:,1:3)));
    diff = max_val - min_val;
    min_val = min_val - diff * 0.1;
    max_val = max_val + diff * 0.1;
    diff = max_val - min_val;
    
    for i = (length(surface)+1) : (length(surface)+1)+ out_num
        surface(i, 1) = min_val + rand * diff;
        surface(i, 2) = min_val + rand * diff;
        surface(i, 3) = min_val + rand * diff;
        surface(i, 4:6) = [255 255 255];
    end
   
    rows = randperm(size(surface,1));
    point_cloud = surface(rows,:);
    %point_cloud = surface;
end

function [inliers, parameters, errors] = RANSAC(points, iterations, err_threshold, num_ind, alg)
    inliers=[];
    parameters = zeros(2, 3);
    indices = 1 : num_ind;
    errors = [];
    
    while(iterations)
        %randomly select + del 
        rows = randperm(size(points,1));
        points = points(rows,:);
        
        %find the parameters -> fit_points()
        switch alg
            case 0 %plane
                [fitted, parameters_tmp] = fit_points_plane(points, indices);
            case 1 %sphere
                [fitted, parameters_tmp] = fit_points_sphere(points, indices);
            case 2 %cilynder
                [fitted, parameters_tmp] = fit_points_cilynder(points, indices);
            otherwise
                inliers = -1;
                parameters = -1;
                errors = -1;
                return;
        end
        
        %
        if(fitted)          
            inliers_tmp = points(1:num_ind, :);
            tmp_ind = num_ind;
            errors_tmp = [];
            
            %count the inlier within a threshold in the remaining set of points
            %-> if fits increase the set -- error under the threshold
            for i = (num_ind+1) : length(points)
                err = error(points(i, 1:3), parameters_tmp);
                errors_tmp(length(errors_tmp) +1) = err;
            end
            
           emin = min(errors_tmp);
           emax = max(errors_tmp);
           lim = emin  + (emax-emin)*err_threshold;
           
            for i = 1 : length(errors_tmp)
                if (errors_tmp(i) < lim) 
                    tmp_ind = tmp_ind + 1;
                    inliers_tmp(tmp_ind, :) = points(i + num_ind, :);
                end
                
            end
            %if more inliers, then earlier -> store it
            if(length(inliers) < length(inliers_tmp))
                inliers = inliers_tmp;
                parameters = parameters_tmp;
                errors = errors_tmp;
            end
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
    mean_p = mean(points);
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
    det = max([det_x det_y det_z]);
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
        fitted = true;
    end
end 

function err = error(point, parameters)
    err = (sum((point - parameters(1,:)) .* parameters(2,:) ))^2;
end

end