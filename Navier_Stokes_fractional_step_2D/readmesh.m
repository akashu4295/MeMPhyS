function [x,y,z] = readmesh(mesh_filename,n)
    opts = delimitedTextImportOptions("NumVariables", 8);
    opts.DataLines = [6, n+5];
    opts.Delimiter = " ";
    opts.VariableNames = ["node_no", "x", "y", "z"];
    opts.SelectedVariableNames = ["node_no", "x", "y", "z"];
    opts.VariableTypes = ["double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";
    
    data = readmatrix(mesh_filename, opts);
    
    x = data(:,2);
    y = data(:,3);
    z = data(:,4);
    clear opts data
end