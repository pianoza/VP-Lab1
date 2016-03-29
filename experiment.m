xl = linspace(-50, 50, 100);      % samples for x;
yl = linspace(-50, 50, 100);      % samples for y
[xx, yy] = meshgrid(xl, yl);    % generating a grid for plotting
z = func(xx, yy);               % getting the values for every sample
figure, mesh(xx, yy, z);        % plotting the levels
axis equal;

rx = Ix;
ry = Iy;
rxx = conv2(rx, dx, 'same');
ryy = conv2(ry, dy, 'same');
rxy = rxx.*ryy;

div = rxx.*ryy-(rxy.^2);
s = (ry.*rxy-rx.*ryy)./div;
t = (rx.*rxy-ry.*rxx)./div;
s(isnan(s)) = 0.0;
t(isnan(t)) = 0.0;
for i = 1:length(co)
    plot(co(i).y+t(co(i).y, co(i).x), co(i).x+s(co(i).y, co(i).x), 'gd', 'MarkerSize',5, 'MarkerFaceColor','k');
end;
hold off;

