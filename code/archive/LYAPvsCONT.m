% clear all; close all; clc;
% 
% I = diag([1 1 1]);
% sym w1 w2 dlt2 x(w1, w2, dlt1) V(x)
% 
% %define state vector x
% x = [w1; w2; dlt1];
% 
% %define the trsnsformation on the state vector x 
% P = I;
% V(x) = x.' * P * x;
% dV = diff(V, x);


%%
I = diag([1 1]);
Q = I;

asd = sin(d);
d_w =..

A = [0,d_w; 3K() 

P = lyap(A,Q);

NP = sqrtm(P);

 = u(NP* J *NP^T)


%%


clear all; close all; clc;

K = 1;
D = 1;
Vg = 1;
E = Vg;
Xs = 1;
Pref = 1;

w_res = 0.1;
d_res = 0.1;
[d,w] = meshgrid(-pi/2:0.1:pi/2, 49:w_res:50);

% eqp = [asin(Pref*Xs/(abs(E)*Vg)), ];
J = [0 1; -3*K*abs(E)*Vg/Xs*cos(d) -K/D];

% x = linspace(-2,2,20);
% y = linspace(-2,2,20);
[x,y] = meshgrid(-5:0.1:5);

f1 = x.^2 + y.^2-10;
f2 = (x+3).^2 + y.^2-10;

hold on;
contourf(x,y,-f1,[0,0]);
contourf(x,y,-f2,[0,0]);
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off;

