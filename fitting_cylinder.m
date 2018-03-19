function fitting_cylinder(points)

cylinder_parameters = cylinder_fitting(points);

end

function cylinder_parameters = cylinder_fitting(points)

        minPC = zeros(3,1);
        minW = zeros(3,1);
        minRSqr = 0;
        minError = 0;
        tmin = 0,
        tmax = 0;
        
        cylinder_parameters = zeros(3); 
        % (1,:) axis origin 
        % (2,:) axis direction
        % (3,1) radius
        % (3,2) height
        % (3,3) min error 
        num_theta_samples  = 1024;
        num_phi_samples = 528;
        num_points = length(points); 
        if (num_points < 6)
            return;
        end

        avg_points = mean(points); 
        centered_points = bsxfun(@minus, points, avg_points);

        [minError, minPC, minW, minRSqr] = fit_by_hemisphere(num_theta_samples, num_phi_samples);

        cylinder_parameters(1,:) = minPC + avg_points;
        cylinder_parameters(2,:) = minPW;
        cylinder_parameters(3,1) = sqrt(minRSqr);
        cylinder_parameters(3,3) = minError;

        for i = 1: num_points
            t = dot(cylinder_parameters(2,:), points(i,:) - cylinder_parameters(1,:));
            tmin = min(t, tmin);
            tmax = max(t, tmax);
        end

        cylinder_parameters(1,:) = cylinder_parameters(1,:) + ( (tmin + tmax) * 0.5 * cylinder_parameters(2,:) );
        cylinder_parameters(3,2) = tmax - tmin;
    
    function [min_error, minPC, minW, minRSqr] = fit_by_hemisphere(num_theta_samples, num_phi_samples)
        iMultiplier = (2*pi) / num_theta_samples;
        jMultiplier = (pi/2) / num_phi_samples;
        minPC= 0;
        minW= [0 0 1];
        minRSqr= 0;
        [min_error, minPC, minRSqr] = error(minW, minPC, minRSqr);
       
        for j = 1 : num_phi_samples
            phi = jMultiplier * j;
            cos_phi = cos(phi);
            sin_phi = sin(phi);
            for i = 1 : num_theta_samples
                theta = iMultiplier * i;
                cos_theta = cos(theta);
                sin_theta = sin(theta);
                
                W = [sin_theta * sin_phi, sin_theta * sin_phi, cos_phi]; %gyanús....
                
                [error_act, PC, rsqr] = error(W);
                if(error_act < min_error)
                    min_error = act_error;
                    minRSqr = rsqr;
                    minW = W;
                    minPC = PC;
                end
            end
        end
    end

    function [error_act, PC, rsqr] = error(W)
        P = eye(3); % - outer_products(W, W); !!
        s = [ 0, -W(2), W(1);
            W(2), 0 , -W(0);
            -W(1), W(0), 0];
        
        Y = length(centered_points); %3*3as mátrxok vektora
        sqrLength  = length(centered_points); %3 elemû vektorok vektora
        
        A = zeros(3);
        B = zeros(3,1);
        qform = 0;
        
        for i = 0 : num_points
            projection = P * centered_points(i,:);
            Y(i, :, :) = projection;
            sqrLength(i, :) = dot(projection, projection);
            A = A + OuterProduct(projection, projection);%!!!
            B = B + sqrLength(i, :) * projection;
            qform = qform + sqrLength(i,:);
        end
        
        
        A = A / num_points;
        B = B / num_points;
        qform = qform / num_points;
        
        Ahat = -S * A * S;
        PC = (Ahat * B) / Trace<Real>(Ahat * A);
        err = 0;
        rsqr = 0;
        
        for  i = 1: num_points
            term = sqrLength(i, :) - Dot(Y(i,:,:), PC) * 2 - qform;
            err = err +  term * term;
            diff = PC - Y(i, :, :);
            rsqr = rsqr + Dot(diff, diff);
        end
  
        error = error / num_points;
        rsqr = rsqr / num_points;
    end
end
