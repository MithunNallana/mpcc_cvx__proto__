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
velocity_max = 5;
velocity_min = 0;
acceleration_max = 1;
acceleration_min = -1;
omega_max = 2;
omega_min = -2;
alpha_max = 1;
alpha_min = -1;

con_weight = 1;
lag_weight = 1;
fwd_weight = 1;
v_smooth_weight = 1;
w_smooth_weight = 1;

%% MPC bro
do_bro = true;
updated = do_bro;

while(do_bro)
    s_guess = s_k + cumsum(velocity_guess*delta_t);
    [ x_path_guess, y_path_guess, theta_path_guess] = progressPath(s_guess, s_k, coff_arc_x, coff_arc_y, arc_param );
    [ coef_x, coef_y, cons_x, cons_y ] = linearizeModel(initial_x, initial_y, initial_theta, velocity_guess, omega_guess, delta_t, n_steps);
    [ Ad,Bd,Cd ] = matrixCostConLag(coef_x, coef_y, cons_x, cons_y, x_path_guess, y_path_guess,theta_path_guess, con_weight, lag_weight, n_steps);
    [ A,B,C ] = matrixCostSmooth( Ad,Bd,Cd,v_smooth_weight, w_smooth_weight, initial_velocity, initial_omega, delta_t, n_steps );
    cvx_begin quiet
    variables v_and_w(2*n_steps,1)
    minimize( v_and_w'*A*v_and_w + B'*v_and_w + C -fwd_weight*sum(v_and_w(1:n_steps,1)*delta_t))
    subject to
    
    v_and_w(1:n_steps,1) >= velocity_min;
    v_and_w(1:n_steps,1) <= velocity_max;
    v_and_w(2:n_steps,1) - v_and_w(1:n_steps-1,1) >= acceleration_min*(delta_t);
    v_and_w(2:n_steps,1) - v_and_w(1:n_steps-1,1) <= acceleration_max*(delta_t);
    v_and_w(1) - initial_velocity >= acceleration_min*(delta_t);
    v_and_w(1) - initial_velocity <= acceleration_max*(delta_t);
%     v_and_w(1:n_steps,1) >= velocity_guess - 0.1;
%     v_and_w(1:n_steps,1) <= velocity_guess + 0.1;
    
    v_and_w(n_steps+1:end,1) >= omega_min;
    v_and_w(n_steps+1:end,1) <= omega_max;
    v_and_w(n_steps+2:end,1) - v_and_w(n_steps+1:end-1,1) >= alpha_min*(delta_t);
    v_and_w(n_steps+2:end,1) - v_and_w(n_steps+1:end-1,1) <= alpha_max*(delta_t);
    v_and_w(n_steps-1) - initial_omega >= alpha_min*(delta_t);
    v_and_w(n_steps-1) - initial_omega <= alpha_max*(delta_t);
%     v_and_w(n_steps+1:end,1) >= omega_guess - 0.05;
%     v_and_w(n_steps+1:end,1) <= omega_guess + 0.05;
    
    cvx_end
    [ x_guess, y_guess, theta_guess ] = calcGuessStates( initial_x, initial_y, initial_theta, velocity_guess, omega_guess, delta_t );
    figure(1);
    axis equal;
    hold on;
    cla;
    plot(x_arc_points(40:60),y_arc_points(40:60),'b');
    plot(x_guess,y_guess,'r');
    plot(x_path_guess,y_path_guess,'g');
    plot(initial_x,initial_y,'bO');
    pause(2);
    velocity_guess = v_and_w(1:n_steps,1);
    omega_guess = v_and_w(n_steps+1:end,1);
    disp(cvx_optval);
    % repeat the loo using updated guess values
%     do_bro = false;
    
end

%% debug
% figure(1);
% hold on;
% plot(x_k,y_k, 'r*');
% plot(initial_x, initial_y, 'b*');