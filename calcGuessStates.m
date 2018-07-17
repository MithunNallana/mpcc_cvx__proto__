function [ x_guess, y_guess, theta_guess ] = calcGuessStates( initial_x, initial_y, initial_theta, velocity_guess, omega_guess, delta_t )
%CALCGUESSSTATES Summary of this function goes here
%   Detailed explanation goes here
theta_diff = omega_guess*delta_t;
theta_guess = cumsum([initial_theta;theta_diff]);
theta_guess = theta_guess(2:end);

vel_diff_cos = velocity_guess .* cos(theta_guess) * delta_t;
vel_diff_sin = velocity_guess .* sin(theta_guess) * delta_t;

x_guess = cumsum([initial_x;vel_diff_cos]);
x_guess = x_guess(2:end);
y_guess = cumsum([initial_y;vel_diff_sin]);
y_guess = y_guess(2:end);

end

