function [ s_k, x_k, y_k  ] = projectToTraj( initial_x, initial_y, x_even,y_even, coff_arc_x, coff_arc_y, arc_param )
%PROJECTTOTRAJ Summary of this function goes here
%   Detailed explanation goes here

cur_index = 1;
cur_dist = Inf;
for i = 1:size(x_even,1)
    if ((x_even(i,1) - initial_x)^2 + (y_even(i,1) - initial_y)^2) < cur_dist
        cur_dist = (x_even(i,1) - initial_x)^2 + (y_even(i,1) - initial_y)^2;
        cur_index = i;
    end
end
cur_index = max(1,cur_index-1);
cur_dist_all = Inf;
for i = cur_index:cur_index+1
   cur_arc_left = arc_param(i);
   x_fun = @(t) (coff_arc_x(i,1)*(t-cur_arc_left)^3 + ...
                 coff_arc_x(i,2)*(t-cur_arc_left)^2 + ...
                 coff_arc_x(i,3)*(t-cur_arc_left) + ...
                 coff_arc_x(i,4));
   y_fun = @(t) (coff_arc_y(i,1)*(t-cur_arc_left)^3 + ...
                 coff_arc_y(i,2)*(t-cur_arc_left)^2 + ...
                 coff_arc_y(i,3)*(t-cur_arc_left) + ...
                 coff_arc_y(i,4));
   num_points = 1000;
   s_cols = linspace(cur_arc_left,(arc_param(2)-arc_param(1)) + cur_arc_left, num_points)';
   for j=1:size(s_cols,1)
       cur_dist = (initial_x - feval(x_fun,s_cols(j)))^2 + ...
                  (initial_y - feval(y_fun,s_cols(j)))^2;
       if cur_dist < cur_dist_all
           cur_dist_all = cur_dist;
           s_k = s_cols(j);
           x_k = feval(x_fun,s_cols(j));
           y_k = feval(y_fun,s_cols(j));
       end
   end
end
end

