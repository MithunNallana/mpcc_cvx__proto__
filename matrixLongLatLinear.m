function [ Ad,Bd,Cd ] = matrixLongLatLinear( initial_x, initial_y, initial_theta, s_k, velocity_guess, omega_guess, path_speed_guess, coff_arc_x, coff_arc_y, arc_param, delta_t, con_weight, lag_weight, n_steps )
%MATRIXLONGLATLINEAR Summary of this function goes here
%   Detailed explanation goes here

[ x_guess, y_guess, theta_guess ] = calcGuessStates( initial_x, initial_y, initial_theta, velocity_guess, omega_guess, delta_t );
s_guess = s_k + cumsum(path_speed_guess*delta_t);
[ x_path_guess, y_path_guess, theta_path_guess, index_guess] = progressPath(s_guess, s_k, coff_arc_x, coff_arc_y, arc_param );

const_con =  (sin(theta_path_guess).*(x_guess-x_path_guess) - cos(theta_path_guess).*(y_guess-y_path_guess))';
const_lag = (-cos(theta_path_guess).*(x_guess-x_path_guess) - sin(theta_path_guess).*(y_guess-y_path_guess))';

% xk sin(o(p))

coef_x = [];
coef_y = [];
utmatrix = triu(ones(n_steps,n_steps));

for m = 1:n_steps
    temp_cos = cos(theta_guess(1:m,1))*delta_t;
    temp_sin = sin(theta_guess(1:m,1))*delta_t;
    x_v_coef = [temp_cos;zeros(n_steps-m,1)];
    y_v_coef = [temp_sin;zeros(n_steps-m,1)];
    temp_u = utmatrix(1:m,1:m);
    x_o_coef = [temp_u*(-temp_sin)*delta_t;zeros(n_steps-m,1)];
    y_o_coef = [temp_u*(+temp_cos)*delta_t;zeros(n_steps-m,1)];
    coef_x = [ coef_x,[ x_v_coef; x_o_coef ] ];
    coef_y = [ coef_y,[ y_v_coef; y_o_coef ] ];
end

coef_x_p = [];
coef_y_p = [];
sin_p = [];
cos_p = [];

for m = 1:n_steps
    temp_index = index_guess(m,1);
    s_pick = arc_param(temp_index,1);
    del_x = 3*coff_arc_x(temp_index,1)*(s_guess(m,1)-s_pick)^2 + ...
            2*coff_arc_x(temp_index,2)*(s_guess(m,1)-s_pick) + ...
            coff_arc_x(temp_index,3);
    del_y = 3*coff_arc_y(temp_index,1)*(s_guess(m,1)-s_pick)^2 + ...
            2*coff_arc_y(temp_index,2)*(s_guess(m,1)-s_pick) + ...
            coff_arc_y(temp_index,3);
    del_del_x = 6*coff_arc_x(temp_index,1)*(s_guess(m,1)-s_pick) + ...
                2*coff_arc_x(temp_index,2);
    del_del_y = 6*coff_arc_y(temp_index,1)*(s_guess(m,1)-s_pick) + ...
                2*coff_arc_y(temp_index,2);
    coef_x_p = [coef_x_p,[del_x*delta_t*ones(m,1);zeros(n_steps-m,1)]];
    coef_y_p = [coef_y_p,[del_y*delta_t*ones(m,1);zeros(n_steps-m,1)]];
    com_diff = (1/(1+(del_y/del_x)^2))*((del_y*del_del_x*delta_t) - (del_x*del_del_y*delta_t))/(del_y^2);
    sin_val = cos(theta_path_guess(m,1))*com_diff;
    cos_val = -sin(theta_path_guess(m,1))*com_diff;
    sin_p = [sin_p,[sin_val*ones(m,1);zeros(n_steps-m,1)]];
    cos_p = [cos_p,[cos_val*ones(m,1);zeros(n_steps-m,1)]];
end

coef_con_new = [];
coef_lag_new = [];

guess_mat = [velocity_guess;omega_guess;path_speed_guess];
for m = 1:n_steps
    coef_con_temp = [sin(theta_path_guess(m,1))*coef_x(:,m);x_guess(m,1)*sin_p(:,m)] - ...
                    [zeros(2*n_steps,1);(sin(theta_path_guess(m,1))*coef_x_p(:,m) + x_path_guess(m,1)*sin_p(:,m))] - ...
                    [cos(theta_path_guess(m,1))*coef_y(:,m);y_guess(m,1)*cos_p(:,m)] + ...
                    [zeros(2*n_steps,1);(cos(theta_path_guess(m,1))*coef_y_p(:,m) + y_path_guess(m,1)*cos_p(:,m))];
    coef_con_new = [coef_con_new,coef_con_temp];
    const_con(1,m) = const_con(1,m) - coef_con_temp'*guess_mat;
    
    coef_lag_temp = -1*[cos(theta_path_guess(m,1))*coef_x(:,m);x_guess(m,1)*cos_p(:,m)] + ...
                    [zeros(2*n_steps,1);(cos(theta_path_guess(m,1))*coef_x_p(:,m) + x_path_guess(m,1)*cos_p(:,m))] - ...
                    [sin(theta_path_guess(m,1))*coef_y(:,m);y_guess(m,1)*sin_p(:,m)] + ...
                    [zeros(2*n_steps,1);(sin(theta_path_guess(m,1))*coef_y_p(:,m) + y_path_guess(m,1)*sin_p(:,m))];
    coef_lag_new = [coef_lag_new,coef_lag_temp];
    const_lag(1,m) = const_lag(1,m) - coef_lag_temp'*guess_mat;
    
end


Ad = zeros(3*n_steps,3*n_steps);
Bd = zeros(3*n_steps,1);
Cd = 0;

for m = 1:n_steps
    Ad = Ad + con_weight*(coef_con_new(:,m)*coef_con_new(:,m)') + lag_weight*(coef_lag_new(:,m)*coef_lag_new(:,m)')  ;
    Bd = Bd + 2*con_weight*coef_con_new(:,m)*const_con(1,m) + 2*lag_weight*coef_lag_new(:,m)*const_lag(1,m);
    Cd = Cd + con_weight*const_con(1,m)*const_con(1,m) + lag_weight*const_lag(1,m)*const_lag(1,m);
end

end

