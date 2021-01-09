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
  CreateProjectionMatrices;
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
function p = CreateProjectionMatrices
  #TODO - how big should the identity matrix be?
  Pn = eye(4);
  Pn = [Pn, zeros(rows(Pn),1)]; 
endfunction

    


