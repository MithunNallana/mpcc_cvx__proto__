function [ x_path_guess, y_path_guess, theta_path_guess, index_guess] = progressPath(s_guess, s_k, coff_arc_x, coff_arc_y, arc_param)
%PROGRESSPATH Summary of this function goes here
%   Detailed explanation goes here
temp_index = 1;
while ( s_k > arc_param(temp_index))
    temp_index = temp_index + 1;
end
temp_index = temp_index - 1;
x_path_guess = zeros(size(s_guess,1),1);
y_path_guess = zeros(size(s_guess,1),1);
index_guess = -1*ones(size(s_guess,1),1);
theta_path_guess = zeros(size(s_guess,1),1);

for i = 1:size(s_guess,1)
    index_mod = min(temp_index+1, size(arc_param,1));
    if( s_guess(i) > arc_param(index_mod))
        temp_index = index_mod;
    end
    s_pick = arc_param(temp_index,1);
    index_guess(i,1) = temp_index;
    x_path_guess(i,1) = coff_arc_x(temp_index,1)*(s_guess(i,1)-s_pick)^3 + ...
                        coff_arc_x(temp_index,2)*(s_guess(i,1)-s_pick)^2 + ...
                        coff_arc_x(temp_index,3)*(s_guess(i,1)-s_pick) + ...
                        coff_arc_x(temp_index,4);
    y_path_guess(i,1) = coff_arc_y(temp_index,1)*(s_guess(i,1)-s_pick)^3 + ...
                        coff_arc_y(temp_index,2)*(s_guess(i,1)-s_pick)^2 + ...
                        coff_arc_y(temp_index,3)*(s_guess(i,1)-s_pick) + ...
                        coff_arc_y(temp_index,4);
    del_x = 3*coff_arc_x(temp_index,1)*(s_guess(i,1)-s_pick)^2 + ...
            2*coff_arc_x(temp_index,2)*(s_guess(i,1)-s_pick) + ...
            coff_arc_x(temp_index,3);
    del_y = 3*coff_arc_y(temp_index,1)*(s_guess(i,1)-s_pick)^2 + ...
            2*coff_arc_y(temp_index,2)*(s_guess(i,1)-s_pick) + ...
            coff_arc_y(temp_index,3);
    theta_path_guess(i,1) = atan2(del_y,del_x);
end

end

