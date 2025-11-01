%% Simple interpolation at specific points
y_query = linspace(1,2,20);
x_query = zeros(size(y_query));
pressure = interpolate_vtk_data('Solution.vtk', x_query, y_query, 'velocity');

plot(y_query,pressure)




%% Interpolation on a grid
% [X, Y] = meshgrid(linspace(-2, 2, 100), linspace(-2, 2, 100));
% P = interpolate_vtk_data('mesh.vtk', X, Y, 'pressure');
% U = interpolate_vtk_data('mesh.vtk', X, Y, 'u');
% V = interpolate_vtk_data('mesh.vtk', X, Y, 'v');
% 
% % Visualize
% contourf(X, Y, P, 20);
% colorbar;