% plot the boundary of the convergence
function plot_contraction(delta_2, omega_2, MU, dist, eq_point)
    hold on;
    MU = real(MU); %remove numerical complex error
    contourf(delta_2, omega_2, MU, [0, 0]);
%     scatter(eq_point(1), eq_point(2));
    contour(delta_2, omega_2, dist, [1, 1]);
%     surf(delta_2, omega_2, MU);
    axis square;
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    hold off;
end