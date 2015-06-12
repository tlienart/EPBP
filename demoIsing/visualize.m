expname = 'demoIsing';
Ngrid   = num2str(200);
Np 		= num2str(50);
Nc 		= num2str(10);
runn    = num2str(1);
lbpd_grid    = load([expname,'_lbpd_grid_np',Ngrid,'.dat']);
lbpd_beliefs = load([expname,'_lbpd_beliefs_np',Ngrid,'.dat']);
lbpd_origbel = load([expname,'_lbpd_origbel_np',Ngrid,'.dat']);
%
epbp_estbel  = load([expname,'_epbp_est_beliefs_np',Np,'_run',runn,'.dat']);
fepbp_estbel = load([expname,'_fepbp_est_beliefs_np',Np,'_nc',Nc,'_run',runn,'.dat']);
%
Np = num2str(50);
pbp_estbel   = load([expname,'_pbp_est_beliefs_np',Np,'_run',runn,'.dat']);

figure
node = 1;
title('Node 1')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
figure
node = 8;
title('Node 7')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
figure
node = 13;
title('Node 16')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
