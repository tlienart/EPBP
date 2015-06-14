expname = 'demoChain';
Ngrid   = num2str(200);
Np 		= num2str(500);
Nc 		= num2str(10);
runn    = num2str(1);
lbpd_grid    = load([expname,'_lbpd_grid_np',Ngrid,'.dat']);
lbpd_beliefs = load([expname,'_lbpd_beliefs_np',Ngrid,'.dat']);
lbpd_origbel = load([expname,'_lbpd_origbel_np',Ngrid,'.dat']);

orig_values = load([expname,'_orig_values.dat']);
obs_values  = load([expname,'_obs_values.dat']);

%
try
    epbp_estbel  = load([expname,'_epbp_est_beliefs_np',Np,'_r',runn,'.dat']);
%    fepbp_estbel = load([expname,'_fepbp_est_beliefs_np',Np,'_nc',Nc,'_r',runn,'.dat']);
catch
    epbp_estbel  = load([expname,'_epbp_est_beliefs_np',Np,'_rrun.dat']);
%    fepbp_estbel = load([expname,'_fepbp_est_beliefs_np',Np,'_nc',Nc,'_rrun.dat']);
end

%
%pbp_estbel   = load([expname,'_pbp_est_beliefs_np',Np,'_r',runn,'.dat']);

figure
node = 1;
title('Node 1')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(orig_values(node)*ones(1,2),linspace(0,0.5,2),'-k')
plot(obs_values(node)*ones(1,2),linspace(0,0.5,2),'-.g')
%plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
%plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
figure
node = 2;
title('Node 2')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(orig_values(node)*ones(1,2),linspace(0,0.5,2),'-k')
plot(obs_values(node)*ones(1,2),linspace(0,0.5,2),'-.g')
%plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
%plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')

figure
node = 3;
title('Node 3')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(orig_values(node)*ones(1,2),linspace(0,0.5,2),'-k')
plot(obs_values(node)*ones(1,2),linspace(0,0.5,2),'-.g')
%plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
%plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
figure
node = 4;
title('Node 4')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(orig_values(node)*ones(1,2),linspace(0,0.5,2),'-k')
plot(obs_values(node)*ones(1,2),linspace(0,0.5,2),'-.g')
%plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
%plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
figure
node = 5;
title('Node 5')
hold on
plot(lbpd_grid,lbpd_beliefs(node,:)/trapz(lbpd_grid,lbpd_beliefs(node,:)),'-b')
plot(lbpd_grid,lbpd_origbel(node,:)/trapz(lbpd_grid,lbpd_origbel(node,:)),'--r')
plot(lbpd_grid,epbp_estbel(node,:)/trapz(lbpd_grid,epbp_estbel(node,:)),'--','color','DarkGreen')
plot(orig_values(node)*ones(1,2),linspace(0,0.5,2),'-k')
plot(obs_values(node)*ones(1,2),linspace(0,0.5,2),'-.g')
%plot(lbpd_grid,fepbp_estbel(node,:)/trapz(lbpd_grid,fepbp_estbel(node,:)),'--','color','Cornflowerblue')
%plot(lbpd_grid,pbp_estbel(node,:)/trapz(lbpd_grid,pbp_estbel(node,:)),'--','color','DarkSalmon')
