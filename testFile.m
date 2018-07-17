%% clearing up stuff
close all;
clear;
clc;

%% Cubic spline generation through waypoints
%
%   In this section of code, we are given a set of waypoints
%   as a input to generate a set of 'natural cubic splines' to
%   between every two waypoints. 'Natural cubic spline' needs two
%   extra initial conditions at initial and final waypoint. 
%   These conditions are (s_first)''(t_0) = 0 and s_last''(t_end) = 0
%   Along with these two constraints we use other continuity constraints
%   of cubicspine at all waypoints to generate curves
%

delta_t = 2;
t_start = 0;
n_points = size(x_points,1);

temp_n = n_points-2;
temp_ones = ones(temp_n,1);
temp_a = diag(4*temp_ones,0);
temp_a = temp_a + diag(temp_ones(1:temp_n-1),-1);
temp_a = temp_a + diag(temp_ones(1:temp_n-1),1);

temp_con = 6/(delta_t)^2;
temp_b_x = temp_con*(x_points(3:end) + x_points(1:end-2) - 2*x_points(2:end-1));
temp_b_y = temp_con*(y_points(3:end) + y_points(1:end-2) - 2*y_points(2:end-1));

temp_sg_x = temp_a\temp_b_x;
temp_sg_y = temp_a\temp_b_y;

temp_sg_x = [0;temp_sg_x;0];
temp_sg_y = [0;temp_sg_y;0];

d_coff_x = x_points(1:end-1);
b_coff_x = temp_sg_x(1:end-1)/2;
a_coff_x = (temp_sg_x(2:end) - temp_sg_x(1:end-1))/(6*delta_t);
c_coff_x = (x_points(2:end)-x_points(1:end-1))/(delta_t);
c_coff_x = c_coff_x - ((delta_t/6)*(2*temp_sg_x(1:end-1)+temp_sg_x(2:end)));
coff_splines_x =  [a_coff_x,b_coff_x,c_coff_x,d_coff_x];

d_coff_y = y_points(1:end-1);
b_coff_y = temp_sg_y(1:end-1)/2;
a_coff_y = (temp_sg_y(2:end) - temp_sg_y(1:end-1))/(6*delta_t);
c_coff_y = (y_points(2:end)-y_points(1:end-1))/(delta_t);
c_coff_y = c_coff_y - ((delta_t*(2*(temp_sg_y(1:end-1))+temp_sg_y(2:end)))/6);
coff_splines_y =  [a_coff_y,b_coff_y,c_coff_y,d_coff_y];

x_spline_points = [];
y_spline_points = [];
segment_lengths = [];
num_points = 100;
length_curve = 0;
segment_lengths = [length_curve];

for i = 1:size(coff_splines_x,1)
    temp_t = linspace(t_start,t_start + delta_t, num_points)';
    for j = 1:size(temp_t,1)
        cur_t = temp_t(j)-t_start;
        cur_x = coff_splines_x(i,1)*cur_t^3 + coff_splines_x(i,2)*cur_t^2 + ...
                coff_splines_x(i,3)*cur_t + coff_splines_x(i,4);
        cur_y = coff_splines_y(i,1)*cur_t^3 + coff_splines_y(i,2)*cur_t^2 + ...
                coff_splines_y(i,3)*cur_t + coff_splines_y(i,4);
        x_spline_points = [x_spline_points;cur_x];
        y_spline_points = [y_spline_points;cur_y];
    end
    s_fun = @(t) sqrt(power((3*coff_splines_x(i,1)*power((t-t_start),2) + ...
                  2*coff_splines_x(i,2)*(t-t_start) + ...
                  coff_splines_x(i,3)),2) + ...
                  power((3*coff_splines_y(i,1)*power((t-t_start),2) + ...
                  2*coff_splines_y(i,2)*(t-t_start) + ...
                  coff_splines_y(i,3)),2));
    length_curve = length_curve + integral(s_fun,t_start,t_start+delta_t);
    segment_lengths = [segment_lengths ; length_curve];
    t_start = t_start + delta_t;
end    


%% arc-length parametrized cubic spline
%   
%   In this section of code, we split parameters of original cubic
%   spline generated earlier in "no_of_ac_parm" points using bisection
%   method such that arc length between these parameters lies the same.
%   We further caluculate (x,y) at these particular points and use 
%   them to generate another set of cubic splines which are now arc
%   length parametrized

no_of_ac_parm = 20;
section_length = length_curve/n_points;
arc_param = linspace(0,length_curve,no_of_ac_parm)';

length_tol = 0.0001;
temp_index = 2;
adjusted_parm = [0];
t_start = 0;
x_even = [x_points(1)];
y_even = [y_points(1)];

for i = 2:no_of_ac_parm-1
    % finding correct spline segment
    while (arc_param(i) >= segment_lengths(temp_index))
        temp_index = temp_index + 1;
    end
    temp_index = temp_index - 1;
    % bisection to do parameter;
    req_length = arc_param(i) - segment_lengths(temp_index);
    t_left = t_start + (temp_index-1)*delta_t;
    t_pick = t_left;
    t_right = t_start + (temp_index)*delta_t;
    t_mid = (t_left + t_right)/2;
    s_fun = @(t) sqrt(power((3*coff_splines_x(temp_index,1)*power((t-t_pick),2) + ...
                  2*coff_splines_x(temp_index,2)*(t-t_pick) + ...
                  coff_splines_x(temp_index,3)),2) + ...
                  power((3*coff_splines_y(temp_index,1)*power((t-t_pick),2) + ...
                  2*coff_splines_y(temp_index,2)*(t-t_pick) + ...
                  coff_splines_y(temp_index,3)),2));
    temp_pss = integral(s_fun,t_pick,t_mid);
    while(abs(temp_pss - req_length) >= length_tol)
        if((temp_pss - req_length) > length_tol)
            t_right = t_mid;
            t_mid = (t_left + t_right)/2;
        else
            t_left = t_mid;
            t_mid = (t_left + t_right)/2;
        end
        temp_pss = integral(s_fun,t_pick,t_mid);
    end
    k_x = (coff_splines_x(temp_index,1)*(t_mid-t_pick)^3 + ...
             coff_splines_x(temp_index,2)*(t_mid-t_pick)^2 + ...
             coff_splines_x(temp_index,3)*(t_mid-t_pick) + ...
             coff_splines_x(temp_index,4));
    k_y = (coff_splines_y(temp_index,1)*(t_mid-t_pick)^3 + ...
             coff_splines_y(temp_index,2)*(t_mid-t_pick)^2 + ...
             coff_splines_y(temp_index,3)*(t_mid-t_pick) + ...
             coff_splines_y(temp_index,4));
    adjusted_parm = [adjusted_parm;t_mid];
    x_even = [x_even;k_x];
    y_even = [y_even;k_y];
end
adjusted_parm = [adjusted_parm; t_start+(n_points-1)*(delta_t)];
x_even = [x_even;x_points(end)];
y_even = [y_even;y_points(end)];

% generation of spline different parameter differences

temp_h = diff(arc_param);
temp_h = temp_h(1);

temp_ones = ones(no_of_ac_parm-2,1);
temp_a = diag(4*temp_ones,0);
temp_a = temp_a + diag(temp_ones(1:no_of_ac_parm-3),-1);
temp_a = temp_a + diag(temp_ones(1:no_of_ac_parm-3),1);

temp_con = 6/(temp_h)^2;
temp_b_x = temp_con*(x_even(3:end) + x_even(1:end-2) - 2*x_even(2:end-1));
temp_b_y = temp_con*(y_even(3:end) + y_even(1:end-2) - 2*y_even(2:end-1));

temp_sg_x = temp_a\temp_b_x;
temp_sg_y = temp_a\temp_b_y;

temp_sg_x = [0;temp_sg_x;0];
temp_sg_y = [0;temp_sg_y;0];

d_coff_x = x_even(1:end-1);
b_coff_x = temp_sg_x(1:end-1)/2;
a_coff_x = (temp_sg_x(2:end) - temp_sg_x(1:end-1))/(6*temp_h);
c_coff_x = (x_even(2:end)-x_even(1:end-1))/(temp_h);
c_coff_x = c_coff_x - ((temp_h/6)*(2*temp_sg_x(1:end-1)+temp_sg_x(2:end)));
coff_arc_x =  [a_coff_x,b_coff_x,c_coff_x,d_coff_x];

d_coff_y = y_even(1:end-1);
b_coff_y = temp_sg_y(1:end-1)/2;
a_coff_y = (temp_sg_y(2:end) - temp_sg_y(1:end-1))/(6*temp_h);
c_coff_y = (y_even(2:end)-y_even(1:end-1))/(temp_h);
c_coff_y = c_coff_y - ((temp_h/6)*(2*(temp_sg_y(1:end-1))+temp_sg_y(2:end)));
coff_arc_y =  [a_coff_y,b_coff_y,c_coff_y,d_coff_y];


%% debugging

x_arc_points = [];
y_arc_points = [];
num_points = 100;
h_start = 0;
temp_hl = diff(arc_param);
temp_hl = temp_hl(1);

for i = 1:size(coff_arc_x,1)
    temp_h = linspace(h_start,h_start + temp_hl, num_points)';
    for j = 1:size(temp_h,1)
        cur_h = temp_h(j)-h_start;
        cur_x = coff_arc_x(i,1)*cur_h^3 + coff_arc_x(i,2)*cur_h^2 + ...
                coff_arc_x(i,3)*cur_h + coff_arc_x(i,4);
        cur_y = coff_arc_y(i,1)*cur_h^3 + coff_arc_y(i,2)*cur_h^2 + ...
                coff_arc_y(i,3)*cur_h + coff_arc_y(i,4);
        x_arc_points = [x_arc_points;cur_x];
        y_arc_points = [y_arc_points;cur_y];
    end
    h_start = h_start + temp_hl;
end


%% Plotting and debugging

figure(1);
hold on;
axis equal;

plot(x_spline_points,y_spline_points,'b');
plot(x_arc_points,y_arc_points,'g');

% plot(x_points,y_points,'r*');
% plot(x_even,y_even,'bO');
