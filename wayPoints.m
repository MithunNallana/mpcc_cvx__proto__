function [ x_points, y_points ] = wayPoints( )
%WAYPOINTGENERATION generates a set of way points from osm maps
%
%   In this version of code way points are given manually
%   Plan is to generate waypoints from osm maps/graph serach
%   based on initial position and goal position

% manual waypoint test data
% x_y_points = [0  , 0  ;
%               100, 0  ;
%               200, 0  ;
%               300, 0  ;
%               400, 200  ;
%               400, 200  ;
%               400, 200  ;
%               400, 200  ;];

% waypoint data from excel
x_y_points = xlsread('x_y_data.xlsx');
x_points = x_y_points(:,1);
y_points = x_y_points(:,2);

end

