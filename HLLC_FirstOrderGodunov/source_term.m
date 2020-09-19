function [var_source_term] = source_term(h, q, tan_theta, cf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% u = q./h;
g = 9.81;
var_source_term = tan_theta*g*h-cf*q.^2./(2*h.^2);

end

