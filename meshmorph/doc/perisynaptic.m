clear all;
close all;
clc;

% suppose synaptic clefts have area a and radius r
% with perisynaptic region of area A and inner and
% outer radii of r and R, respectively
%
% so A_a is ratio of areas A/a
% and
% R_r_minus_1 is ratio of region widths (so to speak)
% = (R-r)/r = R/r-1
 
A_a = [0.1:0.1:10];
R_r_minus_1 = sqrt(1.0+A_a)-1.0;
c = [A_a' R_r_minus_1']


% results:
%   think R_r_minus_1 should be 2
%   therefore A_a = 8
