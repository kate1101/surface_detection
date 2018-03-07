function c = plane_detect()

norm_v = ([0.44 -0.61 0.66]);
scalars = rand(100,3);
outlier_points = rand(50,3);
outlier_points (1:15, 1) = -1 * outlier_points (1:15,1);
outlier_points (10:20, 2)= -1* outlier_points (10:20, 2) ;
outlier_points (25:33, 3) = -1* outlier_points (25:33, 3);

for i=1:100
    plane_points(i,:) = scalars(i, :) .* norm_v;
end
%all_points(i,:) = ([plane_points(1:50,:) ; outlier_points(1:50,:) ; plane_points(51:100,:) ]);
all_points = ([plane_points(1:50,:) ; outlier_points(1:50,:) ; plane_points(51:100,:) ]);
res = pca(all_points);

S = repmat([50,25,10],50,1);
C = repmat([1,2,3],50,1);
s = S(:);
c = C(:);
figure
scatter3(all_points(:,1),all_points(:,2), all_points(:,3), S(:), C(:));
view(40,35)
 a = all_points;
 num = 150;
 r = zeros(num, 1);
g = zeros(num, 1);
b = zeros(num, 1);

c = [a, r, g, b];
end
