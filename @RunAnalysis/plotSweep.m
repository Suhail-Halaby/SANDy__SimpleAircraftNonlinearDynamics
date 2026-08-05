function plotSweep(~, Range, Angle, Cf, Cm)
labels = ["C_X","C_Y","C_Z","C_l","C_m","C_n"];
data   = [Cf; Cm];
xlab   = sprintf("\\%s [deg]", Angle);

figure('Name', "Coefficient Sweep", 'Color', 'w');
for i = 1:6
    subplot(2, 3, i);
    plot(Range, data(i,:), '-o', 'LineWidth', 1.2, 'MarkerSize', 4);
    grid on;
    xlabel(xlab);
    ylabel(labels(i));
    title(labels(i));
end
end