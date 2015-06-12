expname = 'demoIsing';
runn = num2str(1);
Np = num2str(50);
Ngrid = num2str(200);
pbp_evals   = load([expname,'_pbp_evals_np',Np,'_run',runn,'.dat']);
pbp_particles = load([expname,'_pbp_particles_np',Np,'_run',runn,'.dat']);
lbpd_grid    = load([expname,'_lbpd_grid_np',Ngrid,'.dat']);
pbp_est_bel = load([expname,'_pbp_est_beliefs_np',Np,'_run',runn,'.dat']);

%stem(pbp_particles(1,:),pbp_evals(1,:))


plot(lbpd_grid,pbp_est_bel(1,:)/trapz(lbpd_grid,pbp_est_bel(1,:)))