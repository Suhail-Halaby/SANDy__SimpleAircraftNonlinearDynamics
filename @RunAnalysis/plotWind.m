function plotWind(~, Range, Angle, Cw, Cm)
xlab  = sprintf("\\%s [deg]", Angle);
names = ["C_D","C_Y","C_L"];
figure('Name', "Wind-axis coefficients", 'Color', 'w');
for i = 1:3
    subplot(2,2,i);
    plot(Range, Cw(i,:), '-o', 'LineWidth', 1.2, 'MarkerSize', 4);
    grid on; xlabel(xlab); ylabel(names(i)); title(names(i));
end
subplot(2,2,4);   % the classic pitching-moment panel
plot(Range, Cm(2,:), '-o', 'LineWidth', 1.2, 'MarkerSize', 4);
grid on; xlabel(xlab); ylabel("C_m"); title("C_m");
end