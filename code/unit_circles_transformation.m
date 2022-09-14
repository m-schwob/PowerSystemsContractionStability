close all;
clear;
clc;


a = 2.5772;
b = 1.0266;
c = 1.0266;
d = 1.5366;
A = [a b; c d]; %most be inversable
assert(((a * d - b * c) ~= 0), "matrix must be inversable");

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
XYlim = [-4, 4];
linesize = 15;
ledgendstyle = {'Location', 'northwest', 'FontSize', 14, 'Interpreter', 'latex'};
titlestyle = {'fontweight', 'bold', 'fontsize', 18, 'Interpreter', 'latex'};

X2 = plotUnitNorm('2');
X2t = A * X2;

subplot(132); hold on; axis equal;
scatter(X2(1, :), X2(2, :), linesize, '.');
scatter(X2t(1, :), X2t(2, :), linesize, '.');
title('$L_2$', titlestyle{:});
legend('$|x|_{L_2}=1$', '$|Tx|_{L_2}=1$', ledgendstyle{:});
set(gca, 'XLim', XYlim, 'YLim', XYlim);

X1 = plotUnitNorm('1');
X1t = A * X1;

subplot (131); hold on; axis equal;
scatter(X1(1, :), X1(2, :), linesize, '.');
scatter(X1t(1, :), X1t(2, :), linesize, '.');
title('$L_1$', titlestyle{:});
legend('$|x|_{L_1}=1$', '$|Tx|_{L_1}=1$', ledgendstyle{:});
set(gca, 'XLim', XYlim, 'YLim', XYlim);

Xinf = plotUnitNorm('inf');
Xinft = A * Xinf;

subplot (133); hold on; axis equal;
scatter(Xinf(1, :), Xinf(2, :), linesize, '.');
scatter(Xinft(1, :), Xinft(2, :), linesize, '.');
title('$L_\infty$', titlestyle{:});
legend('$|x|_{L_\infty}=1$', '$|Tx|_{L_\infty}=1$', ledgendstyle{:});
set(gca, 'XLim', XYlim, 'YLim', XYlim);

pos = get(gcf, 'Position');
set(gcf, 'Position', pos + [-500 -200 1000 100])

print -depsc unit_circles

%%
function Xn = plotUnitNorm(normType)
    X = randn(10000, 2);

    if (strcmp(normType, 'inf'))
        p_norm = max(abs(X), [], 2);
    end

    if (strcmp(normType, '1'))
        p_norm = sum(abs(X), 2);
    end

    if (strcmp(normType, '2'))
        p_norm = sum(abs(X).^2, 2).^(1/2);
    end

    Xn = bsxfun(@times, X, 1 ./ p_norm);
    Xn = Xn';
    %figure(1); scatter(Xn(:, 1), Xn(:, 2), '.');
end
