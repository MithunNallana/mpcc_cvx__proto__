function [ A,B,C ] = matrixCostConLag(coef_x, coef_y, cons_x, cons_y, x_path_guess, y_path_guess,theta_path_guess, con_weight, lag_weight, n_steps)
%MATRIXCOST Summary of this function goes here
%   Detailed explanation goes here

A = zeros(2*n_steps,2*n_steps);
B = zeros(2*n_steps,1);
C = 0;

for m = 1:n_steps
    con_made = sin(theta_path_guess(m,1))*coef_x(:,m) - ...
                cos(theta_path_guess(m,1))*coef_y(:,m);
    lag_made = -cos(theta_path_guess(m,1))*coef_x(:,m) - ...
                sin(theta_path_guess(m,1))*coef_y(:,m);
    con_cons = sin(theta_path_guess(m,1))*(cons_x(1,m)-x_path_guess(m,1)) - ...
                cos(theta_path_guess(m,1))*(cons_y(1,m)-y_path_guess(m,1));
    lag_cons = -cos(theta_path_guess(m,1))*(cons_x(1,m)-x_path_guess(m,1)) - ...
                sin(theta_path_guess(m,1))*(cons_y(1,m)-y_path_guess(m,1));
    A = A + con_weight*(con_made*con_made') + lag_weight*(lag_made*lag_made');
    B = B + 2*con_weight*con_made*con_cons + 2*lag_made*lag_cons;
    C = C + con_weight*con_cons*con_cons + lag_weight*lag_cons*lag_cons;

end
end

