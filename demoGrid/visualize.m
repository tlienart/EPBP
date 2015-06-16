expname = 'demoGrid';
Ngrid   = num2str(200);
Np 		= num2str(100);
Nc 		= num2str(10);
runn    = num2str(1);
lbpd_grid    = load([expname,'_lbpd_grid_np',Ngrid,'.dat']);
lbpd_beliefs = load([expname,'_lbpd_beliefs_np',Ngrid,'.dat']);
lbpd_origbel = load([expname,'_lbpd_origbel_np',Ngrid,'.dat']);
%
%try
epbp_estbel  = load([expname,'_epbp_est_beliefs_np',Np,'_r',runn,'.dat']);
%    fepbp_estbel = load([expname,'_fepbp_est_beliefs_np',Np,'_nc',Nc,'_r',runn,'.dat']);
%catch
% epbp_estbel  = load([expname,'_epbp_est_beliefs_np',Np,'_rrun.dat']);
%    fepbp_estbel = load([expname,'_fepbp_est_beliefs_np',Np,'_nc',Nc,'_rrun.dat']);
%end
epbp_qmoms = load([expname,'_epbp_qmom_np',Np,'_r',runn,'.dat']);
ep_qmoms   = load([expname,'_ep_qmoments.dat']);
%
%pbp_estbel   = load([expname,'_pbp_est_beliefs_np',Np,'_r',runn,'.dat']);

figure
node = 1;
title('Node 1')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(lbpd_grid,normpdf(lbpd_grid,epbp_qmoms(node,1),epbp_qmoms(node,2)),'-','color','DarkGreen')
plot(lbpd_grid,normpdf(lbpd_grid,ep_qmoms(node,1),ep_qmoms(node,2)),'-','color','Magenta')
%plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
%plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
figure
node = 3;
title('Node 3')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(lbpd_grid,normpdf(lbpd_grid,epbp_qmoms(node,1),epbp_qmoms(node,2)),'-','color','DarkGreen')
plot(lbpd_grid,normpdf(lbpd_grid,ep_qmoms(node,1),ep_qmoms(node,2)),'-','color','Magenta')
%plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
%plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
figure
node = 8;
title('Node 8')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(lbpd_grid,normpdf(lbpd_grid,epbp_qmoms(node,1),epbp_qmoms(node,2)),'-','color','DarkGreen')
plot(lbpd_grid,normpdf(lbpd_grid,ep_qmoms(node,1),ep_qmoms(node,2)),'-','color','Magenta')
%plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
%plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
figure
node = 13;
title('Node 16')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(lbpd_grid,normpdf(lbpd_grid,epbp_qmoms(node,1),epbp_qmoms(node,2)),'-','color','DarkGreen')
plot(lbpd_grid,normpdf(lbpd_grid,ep_qmoms(node,1),ep_qmoms(node,2)),'-','color','Magenta')
%plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
%plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
figure
node = 22;
title('Node 22')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(lbpd_grid,normpdf(lbpd_grid,epbp_qmoms(node,1),epbp_qmoms(node,2)),'-','color','DarkGreen')
plot(lbpd_grid,normpdf(lbpd_grid,ep_qmoms(node,1),ep_qmoms(node,2)),'-','color','Magenta')
%plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
%plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
