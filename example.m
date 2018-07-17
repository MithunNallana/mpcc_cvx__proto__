%% clearing all the data and clearing command window
close all;
clear;
clc;

%% initalization of all required parameters and way points

[initial_x, initial_y, initial_theta, initial_velocity, initial_omega, no_of_ac_curves ] = inputParameters();
[x_points, y_points] = wayPoints();

% construct a parametric equation and try to maintain a parameter
debug_param_plot = 1;
[ coff_arc_x, coff_arc_y, arc_param, x_even, y_even, x_arc_points, y_arc_points ] = parametricPath( x_points, y_points, no_of_ac_curves, debug_param_plot );

% fix initial positions of robot and predict robots motion model x and y
[ s_k, x_k, y_k ] = projectToTraj(initial_x, initial_y, x_even,y_even, coff_arc_x, coff_arc_y, arc_param );

n_steps = 50;
delta_t = 0.1;
velocity_guess = 0.1*ones(n_steps,1);
omega_guess = zeros(n_steps,1);
path_speed_guess = 0.05*ones(n_steps,1);

velocity_max = 5;
velocity_min = 0;
path_speed_max = 7;
path_speed_min = 0;
path_length_max = arc_param(end,1);
path_length_min = arc_param(1,1);

acceleration_max = 1;
acceleration_min = -1;
omega_max = 2;
omega_min = -2;
alpha_max = 1;
alpha_min = -1;

con_weight = 5;
lag_weight = 10000;
fwd_weight = 1;
v_smooth_weight = 1;
w_smooth_weight = 1;
p_smooth_weight = 1;

mat_u = triu(ones(n_steps,n_steps));

%% MPC bro
change_cost = true;
updated = change_cost;

while(change_cost)
    initial_s = s_k;
    [ Ad,Bd,Cd ] = matrixLongLatLinear( initial_x, initial_y, initial_theta, s_k, velocity_guess, omega_guess, path_speed_guess, coff_arc_x, coff_arc_y, arc_param, delta_t, con_weight, lag_weight, n_steps );
    [ A,B,C ] = matrixCostSmooth( Ad,Bd,Cd,v_smooth_weight, w_smooth_weight, p_smooth_weight, initial_velocity, initial_omega, initial_s, delta_t, n_steps );
    cvx_begin quiet
    variables v_and_w(3*n_steps,1)
    minimize( v_and_w'*A*v_and_w + B'*v_and_w + C -fwd_weight*sum(v_and_w(2*n_steps+1:end,1)))
    
    subject to
    v_and_w(1:n_steps,1) >= velocity_min;
    v_and_w(1:n_steps,1) <= velocity_max;
    v_and_w(2:n_steps,1) - v_and_w(1:n_steps-1,1) >= acceleration_min*(delta_t);
    v_and_w(2:n_steps,1) - v_and_w(1:n_steps-1,1) <= acceleration_max*(delta_t);
    v_and_w(1) - initial_velocity >= acceleration_min*(delta_t);
    v_and_w(1) - initial_velocity <= acceleration_max*(delta_t);
%     v_and_w(1:n_steps,1) >= velocity_guess - 0.1;
%     v_and_w(1:n_steps,1) <= velocity_guess + 0.1;
    
    v_and_w(n_steps+1:2*n_steps,1) >= omega_min;
    v_and_w(n_steps+1:2*n_steps,1) <= omega_max;
    v_and_w(n_steps+2:2*n_steps,1) - v_and_w(n_steps+1:2*n_steps-1,1) >= alpha_min*(delta_t);
    v_and_w(n_steps+2:2*n_steps,1) - v_and_w(n_steps+1:2*n_steps-1,1) <= alpha_max*(delta_t);
    v_and_w(n_steps-1,1) - initial_omega >= alpha_min*(delta_t);
    v_and_w(n_steps-1,1) - initial_omega <= alpha_max*(delta_t);
%     v_and_w(n_steps+1:end,1) >= omega_guess - 0.05;
%     v_and_w(n_steps+1:end,1) <= omega_guess + 0.05;

    v_and_w(2*n_steps+1:end,1) >= path_speed_min;
    v_and_w(2*n_steps+1:end,1) <= path_speed_max;
    s_k + mat_u*v_and_w(2*n_steps+1:end,1) >= path_length_min;
    s_k + mat_u*v_and_w(2*n_steps+1:end,1) <= path_length_max;
    cvx_end
    [ x_guess, y_guess, theta_guess ] = calcGuessStates( initial_x, initial_y, initial_theta, velocity_guess, omega_guess, delta_t );
    figure(1);
    axis equal;
    hold on;
    cla;
    plot(x_arc_points(40:60),y_arc_points(40:60),'b');
    plot(x_guess,y_guess,'r');
%     plot(x_path_guess,y_path_guess,'g');
    plot(initial_x,initial_y,'bO');
    pause(2);
    velocity_guess = v_and_w(1:n_steps,1);
    omega_guess = v_and_w(n_steps+1:2*n_steps,1);
    path_speed_guess = v_and_w(2*n_steps+1:end,1);
    disp(cvx_optval);
    % repeat the loo using updated guess values
%     do_bro = false;
    
end

%% debug
% figure(1);
% hold on;
% plot(x_k,y_k, 'r*');
% plot(initial_x, initial_y, 'b*');