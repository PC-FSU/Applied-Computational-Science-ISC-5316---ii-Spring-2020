cloud = load('cloud.dat');

x = -20:0.1:20;
y = -20:0.1:20;
[X, Y] = meshgrid(x, y);
Z = 1 + (sin(X)).^2 + (sin(Y)).^2 - 0.1*exp(-X.^2 -Y.^2);

titles = ["i=1000", "i=5000", "i=7500", "i=10000"];

close all;
for i=1:4
    figure; hold on;
    contour(X, Y, Z);
    scatter(cloud(2*i-1, :), cloud(2*i, :), 'ko');
    xlabel('x');
    ylabel('y');
    title(titles(i));
    hold off;
end

 