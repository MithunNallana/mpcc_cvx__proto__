function [ coef_x, coef_y, cons_x, cons_y ] = linearizeModel( initial_x, initial_y, initial_theta, velocity_guess, omega_guess, delta_t, n_steps)
%LINEARIZEMODEL Summary of this function goes here
%   Detailed explanation goes here

[ x_guess, y_guess, theta_guess ] = calcGuessStates( initial_x, initial_y, initial_theta, velocity_guess, omega_guess, delta_t );

coef_x = [];
coef_y = [];
cons_x = [];
cons_y = [];
utmatrix = triu(ones(n_steps,n_steps));

for m = 1:n_steps
    temp_cos = cos(theta_guess(1:m,1))*delta_t;
    temp_sin = sin(theta_guess(1:m,1))*delta_t;
    x_v_coef = [temp_cos;zeros(n_steps-m,1)];
    y_v_coef = [temp_sin;zeros(n_steps-m,1)];
    temp_u = utmatrix(1:m,1:m);
    x_o_coef = [temp_u*(-temp_sin)*delta_t;zeros(n_steps-m,1)];
    y_o_coef = [temp_u*(+temp_cos)*delta_t;zeros(n_steps-m,1)];
    cons_x = [cons_x, (x_guess(m,1) - x_v_coef'*velocity_guess - x_o_coef'*omega_guess)];
    cons_y = [cons_y, (y_guess(m,1) - y_v_coef'*velocity_guess - y_o_coef'*omega_guess)];
    coef_x = [ coef_x,[ x_v_coef; x_o_coef ] ];
    coef_y = [ coef_y,[ y_v_coef; y_o_coef ] ];
end

end

