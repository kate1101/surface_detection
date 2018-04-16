function fitting_cylinder(file_in)

    points = importdata(file_in);
    points = points(:,1:3);
    cylinder_parameters = cylinder_fitting(points);
    disp('cylinder_parameters:');
    disp(cylinder_parameters);
    
    sphere_radius = cylinder_parameters(3,1);
    sphere_translation = getTrans(cylinder_parameters(1,:));
    sphere_rotation = getRot(cylinder_parameters(2,:)); %% (p1 - rot ?)
    nump = 4000;
    u_max = cylinder_parameters(3,2)*2; % 200;
    m = gen_margin(sphere_translation, cylinder_parameters(2,:), u_max, nump);
    c = gen_points(sphere_radius, sphere_translation, u_max, nump, sphere_rotation);
    save_m('cylinder_nearly.xyz', c);
    save_m('cylinder_margin.xyz', m);
end

function m = gen_margin(t, r, u_max, num)
    %u = u_max * rand(num, 1);    
    %v = 2 * pi * rand(num, 1);
    %x = r .* cos(v) + t(1); 
    %y = r .* sin(v) + t(2);
    %z = u + t(3);
    n =  repmat( (u_max * rand(num, 1)), 3);
    m = zeros(num,3);
    for i = 1: num/2
        m(i,:) = n(i,:) .* r + t;
    end
    for i =  num/2 : num
        m(i,:) = n(i,:) .* -r + t;
    end
    %m(3,:) = m(3,:) - (u_max/2);
    %c = rotZ( thetaXYZ(3), rotY( thetaXYZ(2), rotX(thetaXYZ(1), c) ) );
end

function c = gen_points(r, t, u_max, num, thetaXYZ)
    % Generate cylinder
    % u : R
    % v : [0,2pi]
    u = u_max * (rand(num, 1) - 0.5);    
    v = 2 * pi * rand(num, 1);
    x = r .* cos(v); %+ t(1); 
    y = r .* sin(v); %+ t(2);
    z = u; %+ t(3);
    
    c = [x , y, z];
    c = rotZ( thetaXYZ(3), rotY( thetaXYZ(2), rotX(thetaXYZ(1), c) ) );
    c = c + repmat(t, length(c), 1) ;
    %rgb = repmat(255, num, 3);
    %s = [x , y, z, rgb];
    
    %c = [x, y, z];   
end


function cylinder_parameters = cylinder_fitting(points)

        minPC = zeros(3,1);
        minW = zeros(3,1);
        minRSqr = 0;
        minError = 0;
        tmin = 0;
        tmax = 0;
        
        cylinder_parameters = zeros(3); 
        % (1,:) axis origin 
        % (2,:) axis direction
        % (3,1) radius
        % (3,2) height
        % (3,3) min error 
        num_theta_samples  = 10; %1024;
        num_phi_samples = 10; %528;
        num_points = length(points); 
        if (num_points < 6)
            return;
        end

        avg_points = mean(points); 
        centered_points = bsxfun(@minus, points, avg_points);

        [minError, minPC, minW, minRSqr] = fit_by_hemisphere(num_theta_samples, num_phi_samples);

        cylinder_parameters(1,:) = minPC' + avg_points;
        cylinder_parameters(2,:) = minW;
        cylinder_parameters(3,1) = sqrt(minRSqr);
        cylinder_parameters(3,3) = minError;

        for i = 1: num_points
            t = cylinder_parameters(2,:) *( points(i,:) - cylinder_parameters(1,:))';
            tmin = min(t, tmin);
            tmax = max(t, tmax);
        end

        cylinder_parameters(1,:) = cylinder_parameters(1,:)+ ( (tmin + tmax) * 0.5 * cylinder_parameters(2,:) ); % !!
        cylinder_parameters(3,2) = tmax - tmin;
    
    function [min_error, minPC, minW, minRSqr] = fit_by_hemisphere(num_theta_samples, num_phi_samples)
        iMultiplier = (2*pi) / num_theta_samples;
        jMultiplier = (pi/2) / num_phi_samples;
        minPC= 0;
        minW= [0 0 1];
        minRSqr= 0;
        [min_error, minPC, minRSqr] = error(minW);
       
        for j = 1 : num_phi_samples
            disp('num_phi_samples:');
            disp(j);
            phi = jMultiplier * j;
            cos_phi = cos(phi);
            sin_phi = sin(phi);
            for i = 0 : (num_theta_samples-1)
                disp('num_theta_samples:');
                disp(i);
                theta = iMultiplier * i;
                cos_theta = cos(theta);
                sin_theta = sin(theta);
                
                W = [cos_theta * sin_phi, sin_theta * sin_phi, cos_phi]; 
                [error_act, PC, rsqr] = error(W);
                if(error_act < min_error)
                    min_error = error_act;
                    minRSqr = rsqr;
                    minW = W;
                    minPC = PC;
                end
            end
        end
    end

    function [error_act, PC, rsqr] = error(W)
        P = eye(3) - W'*W; % eye(3); % - outer_products(W, W); !!
        S = [ 0, -W(3), W(2);
            W(3), 0 , (-W(1));
            (-W(2)), W(1), 0];
        
        Y = zeros(length(centered_points),3); %3*3as mátrxok vektora
        sqrLength  = zeros(length(centered_points),1); %3 elemû vektorok vektora
        
        A = zeros(3);
        B = zeros(3,1);
        qform = 0;
        
        for i = 1 : num_points
            projection = P * centered_points(i,:)';
            Y(i, :) = projection;
            sqrLength(i) = projection'*projection;
            A = A + (projection*projection');%!!!
            B = B + sqrLength(i) * projection;
            qform = qform + sqrLength(i);
        end
        
        
        A = A / num_points;
        B = B / num_points;
        qform = qform / num_points;
        
        Ahat = -S * A * S;
        T = trace(Ahat * A);
        if T == 0
            PC = [0; 0; 0];
        else
            PC = (Ahat * B) / T;
        end
        error_act = 0;
        rsqr = 0;
        
        for  i = 1: num_points
            term = sqrLength(i) - (Y(i,:) * PC) * 2 - qform; % PC'
            error_act = error_act +  term * term;
            diff = PC' - Y(i, :);
            rsqr = rsqr + dot(diff, diff);
        end
  
        error_act = error_act / num_points;
        rsqr = rsqr / num_points;
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

function points = rotX(theta, points)
    for i = 1: length(points)
        %X = x;
        Y = points(i,2)*cos(theta) - points(i,3)*sin(theta);
        Z = points(i,2)*sin(theta) + points(i,3)*cos(theta);
        points(i,2) = Y;
        points(i,3) = Z;
    end
end

function points = rotY(theta, points)
    for i = 1: length(points)
        X = points(i,1)*cos(theta) + points(i,3)*sin(theta);
        %Y = y;
        Z = points(i,3)*cos(theta) - points(i,1)*sin(theta);
        points(i,1) = X;
        points(i,3) = Z;
    end
end

function points = rotZ(theta, points)
    for i = 1: length(points)
        X = points(i,1)*cos(theta) - points(i,2)*sin(theta);
        Y = points(i,1)*sin(theta) + points(i,2)*cos(theta);
        %Z = z;
        points(i,1) = X;
        points(i,2) = Y;
    end
end

function sphere_rotation = getRot(vec3)
    y = asin(vec3(1));
    x = atan(vec3(2)/ vec3(3));
    %y = atan( (vec3(1)/ vec3(2)) * sin(x) );  
    %x = asin(vec3(2) / cos(y));
    disp(sin(x)*cos(y) - vec3(2)); 
    disp(cos(x)*cos(y) - vec3(3)); 

    %x = acos ([1 0 0] * vec3');
    %y = acos ([0 1 0] * vec3');
    %z = acos ([0 0 1] * vec3');
    sphere_rotation = [-x y 0];
end

function sphere_translation = getTrans(centere3)
    sphere_translation = centere3;
end