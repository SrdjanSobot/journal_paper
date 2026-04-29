% 
% load('Pe_plot_1.mat');
% plot(SNR_dBm(1,:), edge_matrix(1,:),'-o')
% yticks([1 2 3 4 5])
% 
% hold on
% %plot(SNR_dBm(2,:), edge_matrix(2,:),'-o')
% %plot(SNR_dBm(3,:), edge_matrix(3,:),'-o')
% %plot(SNR_dBm(4,:), edge_matrix(4,:),'-o')
% plot(SNR_dBm(:,6), edge_matrix(:,6),'-o')
% %plot(SNR_dBm(6,:), edge_matrix(6,:),'-o')
% plot(SNR_dBm(7,:), edge_matrix(7,:),'-o')
% %plot(SNR_dBm(8,:), edge_matrix(8,:),'-o')
% %plot(SNR_dBm(9,:), edge_matrix(9,:),'-o')
% %plot(SNR_dBm(10,:), edge_matrix(10,:),'-o')
% hold off;clear; clc;
%% P_e plot
figure;
d = 500:50:8000;
load('plot_h0.mat')
plot(d,P_e_plot, 'LineWidth', 1.5);
grid on;
legend('1 packet per frame','2 packets per frame','3 packets per frame', '4 packets per frame', '5 packets per frame')
xlabel('Distance {\itd} [m]')
ylabel('Error Probability (P_e)')
title('h_{UAV} = h_{BS}        R = k/n * m');

%% plot efficiency
figure;
load("plot_h0.mat")
d = 500:50:8000;
plot(d, eff_edge_, 'LineWidth',1.5);
xlim([500 8000]);
hold on;
plot(0:8000-1,ones(1,8000).*eff_local_(1),'--','LineWidth',1)
%plot(0:10000-1,ones(1,10000).*eff_local(2),'--','LineWidth',1)
%plot(0:10000-1,ones(1,10000).*eff_local(3),'--','LineWidth',1)
%plot(0:10000-1,ones(1,10000).*eff_local(4),'--','LineWidth',1)
%plot(0:10000-1,ones(1,10000).*eff_local(5),'--','LineWidth',1)


grid on;
ylabel('Efficiency (Successful packets / Joule')
xlabel('Distance {\itd} [m]')
title('h_{UAV} = h_{BS}        R = k/n * m');
legend('Edge (m = 1 packet)', 'Edge (m = 2 packets)', 'Edge (m = 3 packets)','Edge (m = 4 packets)','Edge (m = 5 packets)','Local processing per packet')

figure;
plot(d,edge_matrix);
%% plot packet distribution
clear;
%d = load("plot_edge_h0.mat","d");
%d = d.d;
d = 500:50:8000;
packet_distribution_h_0 = load("plot_h0.mat","edge_matrix");
packet_distribution_h_0 = packet_distribution_h_0.edge_matrix;

%packet_distribution_h_50 = load("plot_h50.mat","edge_matrix");
%packet_distribution_h_50 = packet_distribution_h_50.edge_matrix;

packet_distribution_h_100 = load("plot_h100.mat","edge_matrix");
packet_distribution_h_100 = packet_distribution_h_100.edge_matrix;

packet_distribution_h_200 = load("plot_h200.mat","edge_matrix");
packet_distribution_h_200 = packet_distribution_h_200.edge_matrix;

packet_distribution_h_neg20 = load("plot_h-20.mat","edge_matrix");
packet_distribution_h_neg20 = packet_distribution_h_neg20.edge_matrix;

figure;
plot(d, packet_distribution_h_neg20,'LineWidth',1.5);
hold on;
grid on;
plot(d, packet_distribution_h_0,'LineWidth',1.5);
%plot(d, packet_distribution_h_50,'LineWidth',1.5);
%plot(d, packet_distribution_h_100,'LineWidth',1.5);
plot(d, packet_distribution_h_200,'LineWidth',1.5);
xlim([500 8000])
%legend('h_{UAV} = h_{BS} - 20','h_{UAV} = h_{BS}', 'h_{UAV} = h_{BS} + 50', 'h_{UAV} = h_{BS} + 100','h_{UAV} = h_{BS} + 200');
legend('h_{UAV} = h_{BS} - 20','h_{UAV} = h_{BS}', 'h_{UAV} = h_{BS} + 100','h_{UAV} = h_{BS} + 200');
xlabel('Distance {\it d} [m]');
ylabel('Number of packets processed at the edge');

