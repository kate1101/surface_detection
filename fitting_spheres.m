function [sphere_center, sphere_radius, eig_val]  = fitting_sphere(file)
    
    sphere_center = [];
    sphere_radius = [];
    eig_val = [];
    iterations = 1000;
    
    %fdir = 'D:\SZTAKI\Diploma\Pontfelho\dzsudzsabyte\';
    %file = strcat(fdir, fname);
    points = importdata(file);
    
    %eig_val = fitting_sphere_quad_10(points);
    [sphere_center, sphere_radius] = fitting_sphere_lsq(points, iterations);
    %[sphere_center, sphere_radius, eig_val] = fitting_sphere_quad_5(points);
    nump = 100;
    s = gen_points(sphere_center, sphere_radius, nump);
    save_m('sphere_nearly.xyz', s);

function s = gen_points(t, r, num)
    % Generate sphere
    % u: [0;2pi]
    % v: [0;pi]

    u = 2 * pi * rand(num, 1);
    v = pi * rand(num, 1);
    x = r .* sin(v) .* cos(u) + t(1); 
    y = r .* sin(v) .* sin(u) + t(2);
    z = r .* cos(v) + t(3);

    %rgb = repmat(255, num, 3);
    %s = [x , y, z, rgb];
    s = [x, y, z];   
end

function save_m(fname, my_matrix)
    fid = fopen(fname,'wt');
    for i = 1:size(my_matrix,1)
        fprintf(fid,'%g ',my_matrix(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end

function [sphere_center, sphere_radius] = fitting_sphere_lsq(points, iterations)
    %%{
    np = 10;
    r_idxs = datasample(1:length(points), np);
	points = points(r_idxs(1:length(r_idxs)),:);
    %}
    
    points = points(:, 1:3);
    mean_p = mean(points(:,:));
    num_p = length(points);
    
    %sphere_center = mean_p;
    
    % center : the section point of the diverter perpendicular plane of the
    % given points
    % ++ a távolabbi pontok legyenek egy 3as tagjai
    
    sphere_center = find_sphere_center(points);

    tStart = tic;
    
    for j = 1 : iterations
        %%{
        l_avg = 0;
        l_abc_avg = [0 0 0];
        sphere_center_tmp = sphere_center;     
        
        for i = 1 : num_p
            p_dist = sphere_center - points(i,1:3);
            l_avg_act = sqrt(dot(p_dist,p_dist));
            l_avg = l_avg + l_avg_act;
            l_abc_avg = l_abc_avg + ( (sphere_center - points(i,:)) / l_avg_act);
        end
        
        sphere_center = mean_p + (l_avg * l_abc_avg) / num_p^2; 
        sphere_radius = l_avg / num_p;
        
        diff_new = sphere_center - sphere_center_tmp;
        %}                     
        
        %{
        tmp_center = sphere_center;
        len_avg = 0;      
        der_len_avg = [0 0 0];
        
        for i = 1 : num_p
            diff = points(i,:) - sphere_center;
            length_p = sqrt(dot(diff,diff));
            if (length_p > 0 )
                len_avg = len_avg + length_p;
                der_len_avg = der_len_avg + (  (sphere_center - points(i,:)) /  length_p); % (a-x_i) / L_i
                %der_len_avg = der_len_avg - (diff/length_p);
            end
        end
    
        len_avg = len_avg / num_p;
        der_len_avg = der_len_avg / num_p;

        sphere_center = mean_p + len_avg * der_len_avg;
        sphere_radius = len_avg;
        
        diff_new = sphere_center - tmp_center;
        %}
        
        %%{
        zero_m = [0 0 0];
        if(diff_new == zero_m)
            break;
        end
        %}
        

    end
    
 
   
end

function sphere_center = find_sphere_center(points)
    
    tmp_centers = zeros(length(points)-3, 3);
    
    tStart = tic;
    
    for i = 1 : ( length(points)-3 )
        mid_points = zeros(3);
        norm_vecs = zeros(3);
        d = zeros(3,1);
        
        mid_points(1,:) = mean(points(i:(i+1),:));
        mid_points(2,:) = mean(points((i+1):(i+2),:));
        mid_points(3,:) = mean(points((i+2):(i+3),:));
        
        norm_vecs(1,:) = points(i,:) - mid_points(1,:);
        norm_vecs(2,:) = points((i+1),:) - mid_points(2,:);
        norm_vecs(3,:) = points((i+2),:) - mid_points(3,:);
        
        norm_vecs(1,:) = norm_vecs(1,:) / norm(norm_vecs(1,:));
        norm_vecs(2,:) = norm_vecs(2,:) / norm(norm_vecs(2,:));
        norm_vecs(3,:) = norm_vecs(3,:) / norm(norm_vecs(3,:));
        
        d(1) = norm_vecs(1,:) * mid_points(1,:)';
        d(2) = norm_vecs(2,:) * mid_points(2,:)';
        d(3) = norm_vecs(3,:) * mid_points(3,:)';
        
        %tmp_centers(i, :) = inv(norm_vecs) * d;
        
        c = poly(norm_vecs);
        ajd_n = polyvalm(c(1:end-1),norm_vecs);
        tmp_centers(i, :) = (ajd_n / det(norm_vecs)) * d;
    end
    
    sphere_center = median(tmp_centers); % ha félünk az outlierektõl, de pl közeli pontokat nézünk, tehát az in >> out akor vehetjük a mediánját
    
    tElapsed = toc(tStart);
end

end
