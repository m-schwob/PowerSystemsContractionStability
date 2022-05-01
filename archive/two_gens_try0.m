%clc;clear;close all;
syms k1 k2 D1 D2 sd2 a real
T = [k2 k1 xs-1;k2 k1 1;1 0 0];
T1 = T^-1;
J = [-k1/D1 0 -3*k1*a*sd2; 0 -k2/D2 3*k2*a*sd2; -1 1 0];
A = (T^(-1))*J*T

