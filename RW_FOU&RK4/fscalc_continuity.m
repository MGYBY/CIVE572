function [ fs1 ] = fscalc_continuity( ind1, q, dx )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fs1 = -1/dx*(q(ind1+1)-q(ind1));

end

