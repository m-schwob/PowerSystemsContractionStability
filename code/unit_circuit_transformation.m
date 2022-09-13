close all;  
clear;
clc;


a=2.5772;
b=1.0266;
c=1.0266;
d=1.5366;
A = [a b;c d]; %most be inversable
assert(((a*d-b*c) ~= 0), "matrix must be inversable");

% t1 = linspace(-1,1);
% l1_d = abs(t1)-1;
% l1_u = -abs(t1)+1;
% 
% figure(1);
% hold on;
% plot (t1,l1_d,'k')
% plot (t1,l1_u,'k')
% 
% 
% t2 = linspace(0,2*pi);
% l2 = [cos(t2); sin(t2)];

figure;

X2 = plotUnitNorm('2')'; 
X2t = A * X2;

subplot(132); hold on;
scatter(X2(1, :), X2(2, :), '.');
scatter(X2t(1, :), X2t(2, :), '.');
daspect([1 1 1])


X1 = plotUnitNorm('1')'; 
X1t = A * X1;

subplot (131); hold on;
scatter(X1(1, :), X1(2, :), '.');
scatter(X1t(1, :), X1t(2, :), '.');
daspect([1 1 1])


Xinf = plotUnitNorm('inf')'; 
Xinft = A * Xinf;

subplot (133); hold on;
scatter(Xinf(1, :), Xinf(2, :), '.');
scatter(Xinft(1, :), Xinft(2, :), '.');
daspect([1 1 1])


print -depsc unit_circles

%%
function Xn = plotUnitNorm(normType)
X = randn(10000,2);
if(strcmp(normType,'inf'))
    p_norm = max(abs(X), [], 2);
end
if(strcmp(normType,'1'))
    p_norm = sum(abs(X),2);
end
if(strcmp(normType,'2'))
    p_norm = sum(abs(X).^2,2).^(1/2);
end
Xn = bsxfun(@times, X, 1./p_norm);
%figure(1); scatter(Xn(:, 1), Xn(:, 2), '.');
end


