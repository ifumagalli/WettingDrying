function val=temp(vertices,elements,tmp,cna,cn_vecchia,nov)

[dx_cn_vecchia,dy_cn_vecchia]=pdegrad(vertices,elements,cn_vecchia);
dx_cn_vecchia = pdeprtni(vertices,elements,dx_cn_vecchia);
dy_cn_vecchia = pdeprtni(vertices,elements,dy_cn_vecchia);
% [dx_cna,dy_cna]=pdegrad(vertices,elements,cna);
% dx_cna = pdeprtni(vertices,elements,dx_cna);
% dy_cna = pdeprtni(vertices,elements,dy_cna);
val = zeros(nov,1);
dt=0.04;
% val(frontnodes) = (cna(frontnodes)-cn_vecchia(frontnodes))./(dt*dx_cn_vecchia(frontnodes)+eps);
val(tmp) = -(cna(tmp)-cn_vecchia(tmp))./(dt*dx_cn_vecchia(tmp)+eps);
% val(tmp) = -(cna(tmp)-cn_vecchia(tmp))./(dt*dx_cna(tmp)+eps);
ys = unique(vertices(2,:));
for i=1:length(ys)
%     [~,idx] = min(abs(vertices(2,:)-(ymax+ymin)/2));
%     ymed = vertices(2,idx);
%     idxs_ymed = find(vertices(2,:) == ymed);
    idxs = find(vertices(2,:) == ys(i));
% % %     idx_front = max(intersect(idxs,frontnodes));
    idx_front = max(intersect(idxs,tmp));
% % %     idx_min = max(intersect(idxs,tmp));
%     idx_max = max(intersect(idxs,frontnodes))+1;
% % %     idx_max = max(intersect(idxs,firstdrynodes));
    val(idxs(idxs >= idx_front)) = val(idx_front);
% % %     val(idxs((idxs >= idx_min) & (idxs <= idx_max))) = val(idx_min);
end
% function una_med=temp(vertices,elements,boundaries)
% save_path = 'C:\Users\Ivan\Documents\temp_Ivan\';
% load_path = 'C:\Users\Ivan\Dropbox\Elena&Ivan\AnEDP2\Progetto ANEDP2\saves\Ritter_rectang40x5_N192x24_1t5_100\';
% load(strcat(load_path,'dati'))
% 
% tvec = linspace(tspan(1),tspan(2),tspan(3)+1);
% ymin = min(vertices(2,:));
% idxs_ymin = find(vertices(2,:) == ymin);
% ymax = max(vertices(2,:));
% idxs_ymax = find(vertices(2,:) == ymax);
% [~,idx] = min(abs(vertices(2,:)-(ymax+ymin)/2));
% ymed = vertices(2,idx);
% idxs_ymed = find(vertices(2,:) == ymed);
% 
% t = 5;
% [u_ex,h_ex,~] = exact_Stoker_Ritter(vertices,elements,boundaries,save_path,[t t 0],0,3,0,1,1); %<StRi<
% load(strcat(load_path,num2str(t,'%.3f'),'.mat'))
% 
% pfig=figure(3110); set(gcf,'Visible','on'); %<StRi<
%         plot(vertices(1,idxs_ymin),cna(idxs_ymin),vertices(1,idxs_ymax),cna(idxs_ymax), ... %<StRi<
%              vertices(1,idxs_ymed),cna(idxs_ymed),vertices(1,idxs_ymed),sqrt(4*g*h_ex(idxs_ymed))) %<StRi<
%         legend('y = y_m_i_n','y = y_m_a_x','y = y_m_e_d','exact 1D solution') %<StRi<
%         title(strcat('c_n at t = ',num2str(t), ' - Predictor')); %<StRi<
%         %print(pfig,'-deps',strcat(save_path,'cn_wet_ysection',num2str(t,'%.3f'),'apred','.eps')); %<StRi<
%         print(pfig,'-djpeg',strcat(save_path,'cn_wet_ysection',num2str(t,'%.3f'),'apred','.jpeg')); %<StRi<
%         saveas(pfig,strcat(save_path,'cn_wet_ysection',num2str(t,'%.3f'),'apred','.fig')); %<StRi<
% pfig=figure(3111); set(gcf,'Visible','on'); %<StRi<
%         plot(vertices(1,idxs_ymin),una(idxs_ymin),vertices(1,idxs_ymax),una(idxs_ymax), ... %<StRi<
%              vertices(1,idxs_ymed),una(idxs_ymed),vertices(1,idxs_ymed),u_ex(idxs_ymed)) %<StRi<
%         legend('y = y_m_i_n','y = y_m_a_x','y = y_m_e_d','exact 1D solution') %<StRi<
%         title(strcat('u_n at t = ',num2str(t), ' - Predictor')); %<StRi<
%         %print(pfig,'-deps',strcat(save_path,'un_wet_ysection',num2str(t,'%.3f'),'apred','.eps')); %<StRi<
%         print(pfig,'-djpeg',strcat(save_path,'un_wet_ysection',num2str(t,'%.3f'),'apred','.jpeg')); %<StRi<
%         saveas(pfig,strcat(save_path,'un_wet_ysection',num2str(t,'%.3f'),'apred','.fig')); %<StRi<
% pfig=figure(3112); set(gcf,'Visible','on'); %<StRi<
%         plot(vertices(1,idxs_ymin),vna(idxs_ymin),vertices(1,idxs_ymax),vna(idxs_ymax), ... %<StRi<
%              vertices(1,idxs_ymed),vna(idxs_ymed),vertices(1,idxs_ymed),0) %<StRi<
%         legend('y = y_m_i_n','y = y_m_a_x','y = y_m_e_d','exact 1D solution') %<StRi<
%         title(strcat('v_n at t = ',num2str(t), ' - Predictor')); %<StRi<
%         %print(pfig,'-deps',strcat(save_path,'vn_wet_ysection',num2str(t,'%.3f'),'apred','.eps')); %<StRi<
%         print(pfig,'-djpeg',strcat(save_path,'vn_wet_ysection',num2str(t,'%.3f'),'apred','.jpeg')); %<StRi<
%         saveas(pfig,strcat(save_path,'vn_wet_ysection',num2str(t,'%.3f'),'apred','.fig')); %<StRi<
% 
% N=12;
% fc = length(idxs_ymed)/40
% fNyq = fc/2;
% b=fir1(N,0.2);
% una_med = una(idxs_ymed);
% una_filtered=filter(b,fNyq,una(idxs_ymed));
% 
% semi_amp_movavg = 2;
% una_filtered = una_med;
% for iii=semi_amp_movavg+1:length(una_med)-semi_amp_movavg-1
%     una_filtered(iii) = mean(una_med(iii-semi_amp_movavg:iii+semi_amp_movavg));
% end
% 
% b=fir1(N,0.2,'high');
% una_non_filtered=filter(b,fNyq,una(idxs_ymed));
% pfig=figure(1);
% plot(vertices(1,idxs_ymed),una(idxs_ymed),vertices(1,idxs_ymed),una_filtered)%,vertices(1,idxs_ymed),una_non_filtered,vertices(1,idxs_ymed),una_non_filtered+una_filtered)
% legend('originale','filtrata')%,'"residuo"','orig+res')
% title('una filtrata')
% % 
% % [ar]=pdetrg(rppp,rttt);
% % [u_ex,h_ex,frontvelocity_ex(1)] = exact_Stoker_Ritter(rppp,rttt,reee,save_path,[1 1 0],0,3,0,1);
% % wdn_el=pdeintrp(rppp,rttt, h_ex);
% % Vol(1)=sum(wdn_el.*ar);
% %     
% % for i=2:length(tvec)
% %     load(strcat(save_path,num2str(tvec(i),'%.3f'),'.mat'),'vol_corr');
% %     Vol(i) = vol_corr;
% %     [u_ex,h_ex,frontvelocity_ex(i)] = exact_Stoker_Ritter(rppp,rttt,reee,save_path,[tvec(i) tvec(i) 0],0,3,0,1);
% % end
% % 
% % pfig=figure(1); % set(gcf,'Visible','off');
% % plot(tvec,Vol)
% % xlabel('t [s]'); ylabel('volume [m^3]')
% % %print(pfig,'-deps',strcat(save_path,'volume.eps'));
% % print(pfig,'-djpeg',strcat(save_path,'volume.jpeg'));
% % saveas(pfig,strcat(save_path,'volume.fig'));
% % 
% % load(strcat('C:\Users\Ivan\Documents\temp_Ivan\',num2str(tvec(end),'%.3f'),'.mat'));
% % pfig=figure(2); % set(gcf,'Visible','off');
% % plot(tvec,frontvelocity_pred,'.k',tvec,frontvelocity_corr,tvec,frontvelocity_ex,'--r')
% % xlabel('t [s]'); ylabel('front velocity [m/s]')
% % legend('predictor','corrector','exact','Location','SouthWest')
% % %print(pfig,'-deps',strcat(save_path,'frontvelocity.eps'));
% % print(pfig,'-djpeg',strcat(save_path,'frontvelocity.jpeg'));
% % saveas(pfig,strcat(save_path,'frontvelocity.fig'));
% % 
% % save(strcat(save_path,'Vol.mat'), 'Vol');
end