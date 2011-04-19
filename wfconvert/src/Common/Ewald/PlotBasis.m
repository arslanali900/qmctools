b = load ('LPQHI.dat');
r = LPQHI(:,1) - 1;

h10 = b(:,5);
h11 = b(:,6);
h12 = b(:,7);

p1 = plot (r, h10, r, 5*h11, r, 50*h12);
axis ([-1 1 -1.1 1.1])
l2=legend ('h_{j0}', 'h_{j1}*5',  'h_{j2}*50', 0);
set (p1, 'LineWidth', 2.0);
l1 = xlabel ('r-r_j');
set (l1, 'FontSize', 20);
set (l2, 'FontSize', 20);
set (get(gcf, 'CurrentAxes'), 'FontSize', 18);
