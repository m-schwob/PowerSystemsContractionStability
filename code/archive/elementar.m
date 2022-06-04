clear all; close all; clc;

import CalcTool.*
import DisplayTool.*

syms a b c d e f g h k x y z
syms d1 d2 k1 k2 asd

M = [a b c; d e f; g h k];
S = [a b c; b d e; c e f];
D = diag([x y z]);
I = diag([1 1 1]);
J = [-k1/d1 0 -3*k1*asd; 0 -k2/d2 3*k2*asd; -1 1 0]; % our matrix



T = [1 0 0; 0 1 0; 1 1 1];
MM = CalcTool.simmat(J,T);
display_symbolic_pre({MM,'T^{-1}JT='});


