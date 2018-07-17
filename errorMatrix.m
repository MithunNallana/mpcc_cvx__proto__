function [ countouringError, lagError ] = errorMatrix( angle,xd,yd,xk,yk )
%ERRORMATRIX caluculates countouring error and lag error with respect to
%arc length parameterized reference path.
%
%   Initally it caluculates transformation matrices with respect to
%   local frame located on the path and this in inturn caluculates
%   local x and y distances. These distances are countouring error and
%   lag error respectively.

% rotation matrix of local frame

a_R_b    = [ sin(angle) , cos(angle);
            -cos(angle) , sin(angle)];
        
a_P_borg = [ xd; yd ];

a_T_bInv = [ [a_R_b' ; zeros(1,2)] , [ -(a_R_b'*a_P_borg) ; 1] ]; 

localDistances = a_T_bInv * [ xk ; yk ; 1 ];

countouringError = localDistances(1,1);
lagError         = localDistances(2,1);

end

