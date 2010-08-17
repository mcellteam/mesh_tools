clear all;
close all;
clc;

L = 1.0;
a = 1.0;
l1 = [0:0.01:1.0];
l2 = L-l1;
E = a*(1./l1+1./l2);
plot(l1,E)
