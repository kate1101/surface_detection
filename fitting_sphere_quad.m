function [sphere_center, sphere_radius, eig_val]  = fitting_sphere_quad(file)
    
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

function eig_val = fitting_squad_10(points)

num_p = length(points);
inv_num_p = 1 / num_p;
A = zeros(10);

for i = 1:num_p
	x = points(i,0);
	y = points(i,1);
	z = points(i,2);
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    xy = x*y;
    xz = x*z;
    yz = y*z;
    x3 = x*x2;
    xy2 = x*y2;
    xz2 = x*z2;
    x2y = x*xy;
    x2z = x*xz;
    xyz = x*y*z;
    y3 = y*y2;
    yz2 = y*z2;
    y2z = y*yz;
    z3 = z*z2;
    x4 = x*x3;
    x2y2 = x*xy2;
    x2z2 = x*xz2;
    x3y = x*x2y;
    x3z = x*x2z;
    x2yz = x*xyz;
    y4 = y*y3;
    y2z2 = y*yz2;
    xy3 = x*y3;
    xy2z = x*y2z;
    y3z = y*y2z;
    z4 = z*z3;
    xyz2 = x*yz2;
    xz3 = x*z3;
    yz3 = y*z3;

    A(1, 2) = A(1, 1) + x;
    A(1, 3) = A(1, 3) + y;
    A(1, 4) = A(1, 4) + z;
    A(1, 5) = A(1, 5) + x2;
    A(1, 6) = A(1, 6) + y2;
    A(1, 7) = A(1, 7) + z2;
    A(1, 8) = A(1, 8) + xy;
    A(1, 9) = A(1, 9) + xz;
    A(1, 10) = A(1, 10) + yz;
    A(2, 5) = A(2, 5) + x3;
    A(2, 6) = A(2, 6) + xy2;
    A(2, 7) = A(2, 7) + xz2;
    A(2, 8) = A(2, 8) + x2y;
    A(2, 9) = A(2, 9) + x2z;
    A(2, 10) = A(2, 10) + xyz;
    A(3, 6) = A(3, 6) + y3;
    A(3, 7) = A(3, 7) + yz2;
    A(3, 10) = A(3, 10) + y2z;
    A(4, 7) = A(4, 7) + z3;
    A(5, 5) = A(5, 5) + x4;
    A(5, 6) = A(5, 6) + x2y2;
    A(5, 7) = A(5, 7) + x2z2;
    A(5, 8) = A(5, 8) + x3y;
    A(5, 9) = A(5, 9) + x3z;
    A(5, 10) = A(5, 10) + x2yz;
    A(6, 6) = A(6, 6) + y4;
    A(6, 7) = A(6, 7) + y2z2;
    A(6, 8) = A(6, 8) + xy3;
    A(6, 9) = A(6, 9) + xy2z;
    A(6, 10) = A(6, 10) + y3z;
    A(7, 7) = A(7, 7) + z4;
    A(7, 8) = A(7, 8) + xyz2;
    A(7, 9) = A(7, 9) + xz3;
    A(7, 10) = A(7, 10) + yz3;
    A(10, 10) = A(10, 10) + y2z2;
end

A(1, 1) = num_p;
A(2, 2) = A(1, 5);
A(2, 3) = A(1, 8);
A(2, 4) = A(1, 9);
A(3, 3) = A(1, 6);
A(3, 4) = A(1, 10);
A(3, 5) = A(2, 8);
A(3, 8) = A(2, 6);
A(3, 9) = A(2, 10);
A(4, 4) = A(1, 7);
A(4, 5) = A(2, 9);
A(4, 6) = A(3, 10);
A(4, 8) = A(2, 10);
A(4, 9) = A(2, 7);
A(4, 10) = A(3, 7);
A(8, 8) = A(5, 6);
A(8, 9) = A(5, 10);
A(8, 10) = A(6, 9);
A(9, 9) = A(5, 7);
A(9, 10) = A(7, 8);
A(10, 10) = A(6, 7);

for row = 1:10
    for col = 1:10
        A(row, col) = A(col, row);
    end
end

    
for row = 1:10
    for col = 1:10
        A(row, col) = A(row, col) * inv_num_p;
    end
end

eig_val = eig(A);

end

function [sphere_center, sphere_radius, eig_val] = fitting_sphere_quad_5(points)

    num_p = length(points);
    A = zeros(5);
    mean_p = mean(points(:,:));
    
	for i = 1 : num_p
        x = points(i,1) - mean_p(1);
        y = points(i,2) - mean_p(2);
        z = points(i,3) - mean_p(3);
        x2 = x*x;
        y2 = y*y;
        z2 = z*z;
        xy = x*y;
        xz = x*z;
        yz = y*z;
        r2 = x2 + y2 + z2;
        xr2 = x*r2;
        yr2 = y*r2;
        zr2 = z*r2;
        r4 = r2*r2;
        A(1, 5) = A(1, 5) + r2;
        A(2, 2) = A(2, 2) + x2;
        A(2, 3) = A(2, 3) + xy;
        A(2, 4) = A(2, 4) + xz;
        A(2, 5) = A(2, 5) + xr2;
        A(3, 3) = A(3, 3) + y2;
        A(3, 4) = A(3, 4) + yz;
        A(3, 5) = A(3, 5) + yr2;
        A(4, 4) = A(4, 4) + z2;
        A(4, 5) = A(4, 5) + zr2;
        A(5, 5) = A(5, 5) + r4;
    end
    
    A(1, 1) = num_p; 
    A(1, 2) = 0;
    A(1, 3) = 0;
    A(1, 4) = 0;

    for row = 1 : 5
        for col = 1 : 5
            A(col, row) = A(row, col);
        end
    end
    
    eig_val = eig(A);
    
    %{
    for row = 1 : 5
        for col = 1 : 5
            A(row, col) =  A(row, col) / num_p;
        end
    end
    
    eig_val = eig(A);
    %inv = 1 / eig_val(4);
    
    coefficients = zeros(4,1);
    for row = 1 : 4
        coefficients(row) = eig_val(row) / eig_val(4);
    end
    %}
    
    %sphere_center(1:3) = (-0.5) * coefficients(2:4);
    %sphere_radius = sqrt( abs(dot(sphere_center, sphere_center) - coefficients(1)));
    
    %%{
    sphere_center = [];
    sphere_center(1:3) = mean_p(1:3)' - (eig_val(2:4)/(2*eig_val(5)));
    sphere_radius = sqrt(eig_val(2)^2 + eig_val(3)^2 + eig_val(4)^2 - 4*eig_val(1)*eig_val(5)) / (2*eig_val(5));
    %}
end

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
