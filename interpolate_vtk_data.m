function [interp_data] = interpolate_vtk_data(vtk_filename, query_x, query_y, field_name)
% INTERPOLATE_VTK_DATA Read VTK file and interpolate data at query points
%
% Inputs:
%   vtk_filename - Name of the VTK file to read
%   query_x      - Vector of x-coordinates for interpolation
%   query_y      - Vector of y-coordinates for interpolation
%   field_name   - Name of field to interpolate ('velocity', 'pressure', 'u', 'v', 'w', 'p')
%
% Outputs:
%   interp_data  - Interpolated values at query points
%
% Example:
%   x_query = linspace(0, 1, 50);
%   y_query = linspace(0, 1, 50);
%   pressure = interpolate_vtk_data('mesh.vtk', x_query, y_query, 'pressure');

    % Read VTK file
    [nodes, elements, field_data] = read_vtk_file(vtk_filename);
    
    % Extract coordinates
    x = nodes(:, 1);
    y = nodes(:, 2);
    
    % Get the requested field
    if strcmpi(field_name, 'velocity')
        % Return velocity magnitude
        if isfield(field_data, 'velocity')
            u = field_data.velocity(:, 1);
            v = field_data.velocity(:, 2);
            field_values = sqrt(u.^2 + v.^2);
        else
            error('Velocity field not found in VTK file');
        end
    elseif strcmpi(field_name, 'pressure') || strcmpi(field_name, 'p')
        if isfield(field_data, 'pressure')
            field_values = field_data.pressure;
        else
            error('Pressure field not found in VTK file');
        end
    elseif strcmpi(field_name, 'u')
        if isfield(field_data, 'velocity')
            field_values = field_data.velocity(:, 1);
        else
            error('Velocity field not found in VTK file');
        end
    elseif strcmpi(field_name, 'v')
        if isfield(field_data, 'velocity')
            field_values = field_data.velocity(:, 2);
        else
            error('Velocity field not found in VTK file');
        end
    elseif strcmpi(field_name, 'w')
        if isfield(field_data, 'velocity')
            field_values = field_data.velocity(:, 3);
        else
            error('Velocity field not found in VTK file');
        end
    else
        error('Unknown field name: %s', field_name);
    end
    
    % Create triangulation object for interpolation
    tri = triangulation(elements, x, y);
    
    % Create scattered interpolant
    F = scatteredInterpolant(x, y, field_values, 'natural', 'none');
    
    % Interpolate at query points
    interp_data = F(query_x(:), query_y(:));
    
    % Find points outside the domain and set to NaN
    ti = pointLocation(tri, query_x(:), query_y(:));
    interp_data(isnan(ti)) = NaN;
    
    % Reshape if query points were matrices
    if ~isvector(query_x) || ~isvector(query_y)
        interp_data = reshape(interp_data, size(query_x));
    end
end

function [nodes, elements, field_data] = read_vtk_file(filename)
% READ_VTK_FILE Read ASCII VTK unstructured grid file
%
% Outputs:
%   nodes      - Nx3 matrix of node coordinates [x, y, z]
%   elements   - Mx3 matrix of element connectivity (1-indexed)
%   field_data - Structure containing field data

    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    % Initialize outputs
    nodes = [];
    elements = [];
    field_data = struct();
    
    try
        % Skip header lines
        for i = 1:4
            fgetl(fid);
        end
        
        % Read POINTS
        line = fgetl(fid);
        tokens = strsplit(strtrim(line));
        if ~strcmpi(tokens{1}, 'POINTS')
            error('Expected POINTS section');
        end
        num_points = str2double(tokens{2});
        
        % Read node coordinates
        nodes = zeros(num_points, 3);
        for i = 1:num_points
            line = fgetl(fid);
            values = sscanf(line, '%f');
            nodes(i, :) = values(1:3);
        end
        
        % Read CELLS
        line = fgetl(fid);
        while isempty(strtrim(line))
            line = fgetl(fid);
        end
        tokens = strsplit(strtrim(line));
        if ~strcmpi(tokens{1}, 'CELLS')
            error('Expected CELLS section');
        end
        num_cells = str2double(tokens{2});
        
        % Read element connectivity
        temp_elements = cell(num_cells, 1);
        for i = 1:num_cells
            line = fgetl(fid);
            values = sscanf(line, '%d');
            num_nodes_per_elem = values(1);
            temp_elements{i} = values(2:end)' + 1; % Convert to 1-indexed
        end
        
        % Assume triangular elements (you can modify for mixed elements)
        elements = vertcat(temp_elements{:});
        
        % Read CELL_TYPES (skip)
        line = fgetl(fid);
        while isempty(strtrim(line))
            line = fgetl(fid);
        end
        for i = 1:num_cells
            fgetl(fid);
        end
        
        % Read POINT_DATA
        line = fgetl(fid);
        while ischar(line)
            line = strtrim(line);
            if isempty(line)
                line = fgetl(fid);
                continue;
            end
            
            tokens = strsplit(line);
            
            if strcmpi(tokens{1}, 'VECTORS')
                % Read vector field
                field_name = lower(tokens{2});
                vector_data = zeros(num_points, 3);
                for i = 1:num_points
                    line = fgetl(fid);
                    values = sscanf(line, '%f');
                    vector_data(i, :) = values(1:3);
                end
                field_data.(field_name) = vector_data;
                
            elseif strcmpi(tokens{1}, 'SCALARS')
                % Read scalar field
                field_name = lower(tokens{2});
                fgetl(fid); % Skip LOOKUP_TABLE line
                scalar_data = zeros(num_points, 1);
                for i = 1:num_points
                    line = fgetl(fid);
                    scalar_data(i) = sscanf(line, '%f', 1);
                end
                field_data.(field_name) = scalar_data;
                
            end
            
            line = fgetl(fid);
        end
        
    catch ME
        fclose(fid);
        rethrow(ME);
    end
    
    fclose(fid);
end

% Example usage and visualization
function demo_vtk_interpolation()
    % Example: Create a grid of query points
    [X, Y] = meshgrid(linspace(-2, 2, 100), linspace(-2, 2, 100));
    
    % Interpolate pressure field
    try
        P = interpolate_vtk_data('n_1.vtk', X, Y, 'pressure');
        
        figure('Position', [100, 100, 800, 600]);
        
        % Plot interpolated pressure
        subplot(2, 2, 1);
        contourf(X, Y, P, 20, 'LineStyle', 'none');
        colorbar;
        title('Interpolated Pressure');
        xlabel('x'); ylabel('y');
        axis equal tight;
        
        % Interpolate velocity components
        U = interpolate_vtk_data('n_1.vtk', X, Y, 'u');
        V = interpolate_vtk_data('n_1.vtk', X, Y, 'v');
        
        subplot(2, 2, 2);
        contourf(X, Y, U, 20, 'LineStyle', 'none');
        colorbar;
        title('Interpolated u-velocity');
        xlabel('x'); ylabel('y');
        axis equal tight;
        
        subplot(2, 2, 3);
        contourf(X, Y, V, 20, 'LineStyle', 'none');
        colorbar;
        title('Interpolated v-velocity');
        xlabel('x'); ylabel('y');
        axis equal tight;
        
        % Velocity magnitude
        vel_mag = sqrt(U.^2 + V.^2);
        subplot(2, 2, 4);
        contourf(X, Y, vel_mag, 20, 'LineStyle', 'none');
        colorbar;
        title('Velocity Magnitude');
        xlabel('x'); ylabel('y');
        axis equal tight;
        
        % Add streamlines
        hold on;
        step = 5;
        streamslice(X(1:step:end, 1:step:end), ...
                    Y(1:step:end, 1:step:end), ...
                    U(1:step:end, 1:step:end), ...
                    V(1:step:end, 1:step:end), 2);
        hold off;
        
    catch ME
        fprintf('Error: %s\n', ME.message);
        fprintf('Make sure the VTK file exists and has the correct format.\n');
    end
end