clear all; close all; clc;

syms a b c d e f g h k x y z
syms d1 d2 k1 k2 asd

M = [a b c; d e f; g h k];
S = [a b c; b d e; c e f];
D = diag([x y z]);
I = diag([1 1 1]);
J = [-k1/d1 0 -3*k1*asd; 0 -k2/d2 3*k2*asd; -1 1 0]; % our matrix

%% test

J = [-k1/d1 0 -3*k1*asd; 0 -k2/d2 3*k2*asd; -1 1 0]; % our matrix

T = [1 0 1; 0 1 0; 0 0 1];
T =  T * inv([k2 k1 0; k2 k1 -1; 1 0 0]);
JJ = CalcTool.simmat(J,T);

u1 = CalcTool.matmis(JJ, 'L1');
u2 = CalcTool.matmis(JJ, 'L2');
uinf = CalcTool.matmis(JJ, 'Linf');

%display_symbolic_pre({inv(T),'T^{-1}='},{JJ,'T^{-1}JT='},{u1,'\mu_1=max'},{uinf,'\mu_{\infty}=max'},{u2,'\mu_2=\lambda_{max}'})%,{X,''})
display_symbolic_pre({JJ,'T^{-1}JT='},{u1,'\mu_1=max'},{uinf,'\mu_{\infty}=max'});

%% symetric test k1=k2, d1=d2

JS = [-k/d 0 -3*k*asd; 0 -k/d 3*k*asd; -1 1 0]; % our matrix

T = [1 0 0; 0 1 0; 1 0 0];
JJ = CalcTool.simmat(JS,T);

u1 = CalcTool.matmis(JJ, 'L1');
u2 = CalcTool.matmis(JJ, 'L2');
uinf = CalcTool.matmis(JJ, 'Linf');

%display_symbolic_pre({inv(T),'T^{-1}='},{JJ,'T^{-1}JT='},{u1,'\mu_1=max'},{uinf,'\mu_{\infty}=max'},{u2,'\mu_2=\lambda_{max}'})%,{X,''})
DisplayTool.display_symbolic_pre({JJ,'T^{-1}JT='},{u1,'\mu_1=max'},{uinf,'\mu_{\infty}=max'});

%% jordan half symetric test d1=d2

%J = [-k1/d1 0 -3*k1*asd; 0 -k2/d2 3*k2*asd; -1 1 0]; % our matrix
J = [-k1/d1 0 -3*k1*asd; 0 -k2/d1 3*k2*asd; -1 1 0];

DisplayTool.display_symbolic_pre({jordan(J),''});

%% jordan symetric test k1=k2, d1=d2

%J = [-k1/d1 0 -3*k1*asd; 0 -k2/d2 3*k2*asd; -1 1 0]; % our matrix
JS = [-k/d 0 -3*k*asd; 0 -k/d 3*k*asd; -1 1 0]; 
[P,Jj] = jordan(JS);
DisplayTool.display_symbolic_pre({Jj,'J_j='},{P,'P='});
