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

%% 
Q = (J.')*S + S*J;
vars = [a b c d e f];
eqn = [Q(1,2)==0,Q(1,3)==0,Q(2,3)==0,Q(1,1)==-1,Q(2,2)==-1,Q(3,3)==-1];
s = solve(eqn, vars);
solvec = [s.a s.b s.c s.d s.e s.f];
DisplayTool.display_symbolic_pre({Q,'-Q='},{solvec,strcat(latex(vars),'=')});
DisplayTool.display_symbolic(s.a, s.b, s.c, s.d, s.e, s.f);

%% symetric system
JS = [-k1/d1 0 -3*k1*asd; 0 -k1/d1 3*k1*asd; -1 1 0]; % our matrix
Q = (JS.')*S + S*JS;
vars = [a b c d e f];
eqn = [Q(1,2)==0,Q(1,3)==0,Q(2,3)==0,Q(1,1)==-1,Q(2,2)==-1,Q(3,3)==-1];
s = solve(eqn, vars);
solvec = [s.a s.b s.c s.d s.e s.f];
DisplayTool.display_symbolic_pre({Q,'-Q='},{solvec,strcat(latex(vars),'=')});

%calculate with the result - not working
x = [(-1-6*asd*k1)*d1/k1, d1/k1, 1, d1/k1*(6*asd*k1-1), -1 , -(18*k1*(asd^2)*(d1^2)-6*asd*(d1^2)+k1)/d1];
% T = [s.a s.b s.c; s.b s.d s.e; s.c s.e s.f];
T = [x(1) x(2) x(3); x(2) x(4) x(5); x(3) x(5) x(6)];
JJ = CalcTool.simmat(JS,T);
DisplayTool.display_symbolic_pre({JJ,'JJ='});


u1 = CalcTool.matmis(JJ, 'L1');
u2 = CalcTool.matmis(JJ, 'L2');
uinf = CalcTool.matmis(JJ, 'Linf');
DisplayTool.display_symbolic_pre({JJ,'T^{-1}JT='},{u1,'\mu_1=max'},{uinf,'\mu_{\infty}=max'});

%% additinal
JS = [-k1/d1 0 -3*k1*asd; 0 -k1/d1 3*k1*asd; -1 1 0]; % our matrix

%trying to inverse P
a=-1+6*k1*asd;
b=1;
c=k1/d1;
d=6*k1*asd-1;
e=-k1/d1;
f=-18*(k1^2)*(asd^2)+6*k1*asd-(k1^2)/(d1^2);

T = [a b c; b d e; c e f];
DisplayTool.display_symbolic_pre({T,'T='},{inv(T),'T^{-1}='});

JJ = CalcTool.simmat(JS,T);
DisplayTool.display_symbolic_pre({JJ,'JJ='});





