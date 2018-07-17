function [ A,B,C ] = matrixCostSmooth( A,B,C,v_smooth_weight, w_smooth_weight, p_smooth_weight, initial_velocity, initial_omega, initial_s, delta_t, n_steps )
%MATRIXCOSTSMOOTH Summary of this function goes here
%   Detailed explanation goes here

temp_v = diag(-1*ones(n_steps,1)) + diag(ones(n_steps-1,1),1);
temp_v = temp_v'*temp_v;
temp_v(1,1) = 2;

A = A + [v_smooth_weight*temp_v,zeros(n_steps,2*n_steps); ...
        zeros(n_steps,n_steps),w_smooth_weight*temp_v,zeros(n_steps,n_steps); ...
        zeros(n_steps,2*n_steps),p_smooth_weight*temp_v]/(delta_t)^2;
    
B = B + [-2*v_smooth_weight*initial_velocity;zeros(n_steps-1,1); ...
        -2*w_smooth_weight*initial_omega;zeros(n_steps-1,1); ...
        -2*p_smooth_weight*initial_s;zeros(n_steps-1,1);]/delta_t;
    
C = C + v_smooth_weight*initial_velocity^2/delta_t^2 + ... 
        w_smooth_weight*initial_omega^2/delta_t^2 + ...
        p_smooth_weight*initial_s^2/delta_t^2;
end

