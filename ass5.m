function ass5 ()
  fh = fopen('bh.dat', 'r');
  A = fscanf(fh, '%f%f%f%f', [4 inf]); fclose(fh);
  first_x = A(1:1, :);
  first_y = A(2:2, :);
  second_x = A(3:3, :);
  second_y = A(4:4, :);
  
  # -----------normalizing image points----------------
  tx_first = mean(first_x);
  ty_first = mean(first_y);
  tx_second = mean(second_x);
  ty_second = mean(second_y);
  
  first_x_cond = first_x - tx_first;
  first_y_cond = first_y - ty_first;
  second_x_cond = second_x - tx_second;
  second_y_cond = second_y - ty_second;
  
  sx_first = mean(abs(first_x_cond));
  sy_first = mean(abs(first_y_cond));
  sx_second = mean(abs(second_x_cond));
  sy_second = mean(abs(second_y_cond));
  
  first_x_cond = first_x_cond/sx_first;
  first_y_cond = first_y_cond/sy_first;
  second_x_cond = second_x_cond/sx_second;
  second_y_cond = second_y_cond/sy_second;
  
  first = [first_x_cond; first_y_cond];
  second = [second_x_cond; second_y_cond];
  # --------------Formulating a homogeneous linear equation----------
  a = [];
  for i = 1 : size(first, 2)
    temp = [first(1,i)*second(1,i) first(2,i)*second(1,i) second(1,i) first(1,i)*second(2,i) first(2,i)*second(2,i) second(2,i) first(1,i) first(2,i) 1];
    a = [a;temp];
    end
   # --------------SVD----------
  [u,s,v] = svd(a);
  f = v(:,9);
  f = transpose(reshape(f,3,3));
  
  t_first = CreateTransformationMatrix(tx_first,ty_first,sx_first,sy_first);
  t_second = CreateTransformationMatrix(tx_second,ty_second,sx_second,sy_second);
  
 
  # --------- Enforcing the internal constraint ----------------
  [u,s,v] = svd(f);
  s(3,3) = 0;
    
  f = u * s * transpose(v);
  
  f = transpose(t_second) * f;
  f = f * t_first
  
  d = det(f)
  
  #------ geometric_error ----------
  geometric_error = geom_dist(f, first_x, first_y, second_x, second_y)
  
  
  #---- ass 5------
  #---- part 1------  
  p = CreateProjectionMatrix_N;
  p_prime = CreateProjectionMatrix_P(f,second_x, second_y);
  X = linear_triangulation(p, p_prime, first_x, first_y, second_x, second_y);
  figure(1); 
  scatter3(X(1,:), X(2,:), X(3,:), 10, 'filled'); 
  axis square; 
  view(32, 75);
  
  #---- part 2------
  
  gh = fopen('pp.dat', 'r');
  B = fscanf(gh, '%f%f%f%f%f%f%f', [7 inf]); fclose(gh);
  pp_first_x = B(1:1, :);
  pp_first_y = B(2:2, :);
  pp_second_x = B(3:3, :);
  pp_second_y = B(4:4, :);
  XE = B(5:5, :);
  YE = B(6:6, :);
  ZE = B(7:7, :);
  
  % triangulate projective object points X2 from the five homologous image points using projection matrices p and p_prime.
  X2 = linear_triangulation(p, p_prime, pp_first_x, pp_first_y, pp_second_x, pp_second_y);
  
  % spatial 3D homography H - determination of spatial transformation of projective object points X2 to 
  % corresponding Euclidean object points XE.
  H = homography3(X2, XE);
  
  % Apply transformation H to all object points of projective reconstruction X
  result = [];
  for point_num = 1 : size(X, 2)
    x = X(:, i) * H;
    result = [result, x];
  end
  
  % normalizing result  
  result(1, 1) /= result(4, 1); 
  result(2, 1) /= result(4, 1); 
  result(3, 1) /= result(4, 1); 
  result(4, 1) /= result(4, 1); 
  
  % visualizing result  
  figure(2); 
  scatter3(result(1,:), result(2,:), result(3,:), 10, 'filled'); 
  axis square; 
  view(32, 75);
  
endfunction

function sum = geom_dist(f, first_x, first_y, second_x, second_y)
  sum = 0;
  for point_number = 1:8
    point = [first_x(point_number), first_y(point_number), 1];
    corr_point = [second_x(point_number), second_y(point_number), 1];
    epipolar_line = f * transpose(point);
    epipolar_line_T =  transpose(f) * transpose(corr_point) ;
    gd = dist(epipolar_line, corr_point)^2 + dist(epipolar_line_T,point)^2;
    sum = sum + gd;
  end
endfunction

function d = dist(l,p)
  d = (l(1)*p(1)+l(2)*p(2)+l(3)*p(3))/(p(2)*sqrt(l(1)^2+l(2)^2));
endfunction

function norm_mat = normalize_image_points(x)

# Calculate center of mass
centroid = mean(x,2);

# Calculate distances of each point to the centroid
dist = sqrt(sum((x - centroid) .^ 2));

# Calculate current mean distance
mean_dist = mean(dist);

#shift the image coordinate systems to the respective centroids
#point coodinates in new cos: (x – centoid-x, y – centroid-y)
x = x-centroid;

# scale so that the mean distance from the origin to a point equals sqrt(2)
#multiplying old coordinates by sqrt(2) and dividing out the mean distance
norm_mat = x * sqrt(2) / mean_dist;

endfunction


function t = CreateTransformationMatrix(tx, ty, sx, sy)
    t_scale = eye(3);
  t_scale(1,1) = 1/sx;
  t_scale(2,2) = 1/sy;
  
  t_translate = eye(3);
  t_translate(1,3) = -tx;
  t_translate(2,3) = -ty;
  
  t = t_scale*t_translate;
endfunction

# caluclate projection matrix N for simple camera 1
function Pn = CreateProjectionMatrix_N
  Pn = eye(3);
  Pn = [Pn, zeros(rows(Pn),1)]; 
endfunction


# caluclate projection matrix P' for second camera
function p = CreateProjectionMatrix_P(F, second_x, second_y)
  # calculate second epipole (by finding intersection of two epipolar lines)
  corr_point_1 = [second_x(1), second_y(1), 1];
  corr_point_2 = [second_x(2), second_y(2), 1];
  epipolar_line_T_1 =  transpose(F) * transpose(corr_point_1);
  epipolar_line_T_2 =  transpose(F) * transpose(corr_point_2);
  e = cross(epipolar_line_T_1,epipolar_line_T_2)
  #skew symmetric matrix derived from second epipole
  skew = [0,-e(3), e(2); e(3), 0, -e(1); -e(2), e(1), 0]
  p = skew * F;
  #add second epipole as 4th column
  p = [p,e]
endfunction

function X_mat = linear_triangulation(p, p_prime, first_x, first_y, second_x, second_y)
  X_mat = [];
  #find  4x4 matrix A for each point pair
  for point_number = 1 : size(first_x, 2)
    A = [first_x(point_number) * p(3, :) - p(1, :); 
        first_y(point_number) * p(3, :) - p(2, :); 
        second_x(point_number) * p_prime(3, :) - p_prime(1, :);
        second_y(point_number) * p_prime(3, :) - p_prime(2, :)
        ];
    # find object point X using constraint AX = 0
    [u,s,v] = svd(A);
    X = v(: ,4);
    X(1, 1) /= X(4, 1); 
    X(2, 1) /= X(4, 1); 
    X(3, 1) /= X(4, 1); 
    X(4, 1) /= X(4, 1); 
    X_mat = [X_mat, X];  
  end
endfunction

function H = homography3(x1, x2)
  T1 = condition3(x1); c1 = T1 * x1;
  T2 = condition3(x2); c2 = T2 * x2;
  A = design_homo3(c1, c2);
  h = solve_dlt(A);
  H = inv(T2) * reshape(h, 4, 4)' * T1;
endfunction

function T = condition3(x)
  % error occurs when x = XE because XE only has 1 line 
  %- unfortunatly we didn't manage to fix this problem (therefore no second figure is portrayed)
  tx = mean(x(1, :));             ty = mean(x(2, :));             tz = mean(x(3, :));
  sx = mean(abs(x(1, :) - tx));   sy = mean(abs(x(2, :) - ty));   sz = mean(abs(x(3, :) - tz));
  T = [1/sx 0     0     -tx/sx;
      0     1/sy  0     -ty/sy;
      0     0     1/sz  -tz/sz;
      0     0     0     1];
endfunction
function A = design_homo3(x1, x2)
  A = [];
  for i = 1 : size(x1, 2)
     A = [A;
          -x2(4, i) * x1(:, i)'   0   0   0   0   0   0   0   0   x2(1, i) * x1(:, i)';
          0   0   0   0   -x2(4, i) * x1(:, i)'   0   0   0   0   x2(2, i) * x1(:, i)';
          0   0   0   0   0   0   0   0   -x2(4, i) * x1(:, i)'   x2(3, i) * x1(:, i)';];
  endfor
endfunction

function x = solve_dlt(A)
  [U, D, V] =  svd(A);
  x = V(:, end);
endfunction
