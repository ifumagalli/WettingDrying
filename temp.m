pause off; close all; our_swesimodeladapt_gio_bc_predcorr(rppppp,rttttt,reeeee,tspan,tstart,theta);
% save_path = 'C:\Users\Ivan\Documents\temp_Ivan\';
% tvec = linspace(1,4.4,85+1);
% 
% [ar]=pdetrg(rppp,rttt);
% [u_ex,h_ex,frontvelocity_ex(1)] = exact_Stoker_Ritter(rppp,rttt,reee,save_path,[1 1 0],0,3,0,1);
% wdn_el=pdeintrp(rppp,rttt, h_ex);
% Vol(1)=sum(wdn_el.*ar);
%     
% for i=2:length(tvec)
%     load(strcat(save_path,num2str(tvec(i),'%.3f'),'.mat'),'vol_corr');
%     Vol(i) = vol_corr;
%     [u_ex,h_ex,frontvelocity_ex(i)] = exact_Stoker_Ritter(rppp,rttt,reee,save_path,[tvec(i) tvec(i) 0],0,3,0,1);
% end
% 
% pfig=figure(1); % set(gcf,'Visible','off');
% plot(tvec,Vol)
% xlabel('t [s]'); ylabel('volume [m^3]')
% %print(pfig,'-deps',strcat(save_path,'volume.eps'));
% print(pfig,'-djpeg',strcat(save_path,'volume.jpeg'));
% saveas(pfig,strcat(save_path,'volume.fig'));
% 
% load(strcat('C:\Users\Ivan\Documents\temp_Ivan\',num2str(tvec(end),'%.3f'),'.mat'));
% pfig=figure(2); % set(gcf,'Visible','off');
% plot(tvec,frontvelocity_pred,'.k',tvec,frontvelocity_corr,tvec,frontvelocity_ex,'--r')
% xlabel('t [s]'); ylabel('front velocity [m/s]')
% legend('predictor','corrector','exact','Location','SouthWest')
% %print(pfig,'-deps',strcat(save_path,'frontvelocity.eps'));
% print(pfig,'-djpeg',strcat(save_path,'frontvelocity.jpeg'));
% saveas(pfig,strcat(save_path,'frontvelocity.fig'));
% 
% save(strcat(save_path,'Vol.mat'), 'Vol');