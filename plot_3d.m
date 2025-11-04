data = readmatrix("Solution.csv");
x = data(:,1);
y = data(:,2);
z = data(:,3);

u = scatteredInterpolant(x,y,z,data(:,4));
v = scatteredInterpolant(x,y,z,data(:,5));
[X,Y] = meshgrid(linspace(-0.5,0.5,50),linspace(-0.5,0.5,50));
Z = zeros(size(Y));

U = u(X,Y,Z);
V = v(X,Y,Z);
contourf(X,Y,U)
hold on
streamslice(X,Y,U,V)
