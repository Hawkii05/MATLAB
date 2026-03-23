function Error(harp_estimates)

flight_traj = harp_estimates.data_out.carp.relative_traj;
outputs = harp_estimates.outputs;

% Extract key positions
pi_x = ft2m(outputs.harp.position_x); %inputs.altitude.landing_location(:,1);
pi_y = ft2m(outputs.harp.position_y); %inputs.altitude.landing_location(:,2);
pi_z = 0; % Ground level

% ERROR PLOT: PI vs Flight Trajectory End

figure('Position', [100, 100, 800, 500]);

% Compute error vector
err = [pi_x, pi_y, pi_z] - flight_traj(end, :);
err_labels = {'East (X)', 'North (Y)', 'Altitude (Z)'};

% Bar chart
b = bar(err, 'FaceColor', 'flat');
b.CData(1,:) = [0.2, 0.5, 0.9];   % Blue  - East
b.CData(2,:) = [0.2, 0.75, 0.3];  % Green - North
b.CData(3,:) = [0.9, 0.4, 0.2];   % Red   - Altitude

% Zero reference line
yline(0, 'k--', 'LineWidth', 1.2);

% Annotate each bar with its value
for i = 1:3
    offset = sign(err(i)) * max(abs(err)) * 0.04;
    text(i, err(i) + offset, sprintf('%.2f m', err(i)), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 11, 'FontWeight', 'bold');
end

% Labels & formatting
set(gca, 'XTickLabel', err_labels, 'FontSize', 11);
xlabel('Component', 'FontSize', 12);
ylabel('Error (m)', 'FontSize', 12);
title(sprintf('PI vs Flight Trajectory — Position Error'), ...
    'FontSize', 13, 'FontWeight', 'bold');
grid on;
box on;
end