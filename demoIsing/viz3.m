grid = load('demoIsing_grid.dat');
epbp = load('demoIsing_epbp_est_beliefs_np50_run1.dat');

plot(grid,epbp(1,:)/trapz(grid,epbp(1,:)))