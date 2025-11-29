function [h_uav, SNR_mW, SNR_dBm] = snr_from_excess_PL_novo(L)
% Model from paper: " Al-Hourani, A. and Gomez, K., 2017. Modeling 
% cellular-to-UAV path-loss for suburban environments. IEEE Wireless 
% Communications Letters, 7(1), pp.82-85."
alfa  = 3.9; % Terrestrial path-loss exponent *Different from the original value in paper*
teta_0 = -3.61; % Angle offset
B = 4.14; % Angle scaler
eta_0 = 20.70; % Excess path-loss offest
a = -0.41; % UAV shafowing slope
sigma_0 = 5.86; % UAV shadowing offset
A = -23.29; % Excell path-loss scaler


%% Uncomment this part for: changing distance from BS vs fixed UAV altitude at -15 0 30 60 90
 d = 500:500:8000;
 %d = 2000;
 h_uav = [L - 40];
% h_uav = [0];

 depression_angle = zeros(length(d), length(h_uav));
 for i = 1:length(d)
     for j = 1:length(h_uav)
        depression_angle(i,j) = acosd( sqrt( d(i)^2 /( d(i)^2 + h_uav(j)^2 ) ) );        
     end    
 end
 depression_angle(:,1) = depression_angle(:,1) * -1;

%% Uncomment this part for: Fixed distances from BS {500 1000 1500 2000 2500}, changing UAV altitude
%Extracted from depression_angle table (in previous section)
   % % d = [500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000];
   %  d = [500 1000 1500 2000 3000 4000 5000 6000]
   %  h_uav = linspace(-20,100,30);
   %  depression_angles = [-2.290610042638435 11.309932474020195;... %500
   %                      -1.145762838175144 5.710593137499632;... % 1000
   %                      -0.763898460929882 3.814074834290378;... % 1500
   %                      -0.572938697683265 2.862405226111779;... % 2000m
   %                      %-0.458356458000046 2.290610042638435;... % 2500m
   %                      -0.381966204729653 1.909152432996390;... % 3000m
   %                     % -0.327400890844405 1.636577041616692;... % 3500m
   %                      -0.286476510277181 1.432096184164545;... % 4000m
   %                      %-0.254646232272433 1.273030020056664;... % 4500m
   %                      -0.229181895753578 1.145762838175144;... % 5000m
   %                      %-0.208347370806743 1.041626676010165;... % 5500m
   %                      -0.190985224360512 0.954841253872355]... % 6000m
   %                      %-0.176294149842681 0.881403996582091;... % 6500m
   %                      %-0.163701781733923 0.818455461688598;... % 7000m
   %                      %-0.152788383203811 0.763898460929882;... % 7500m
   %                      % -0.143239150369941 0.716159945470150];   % 8000m
   %  %depression_angles = [  -0.381966204729653 1.909152432996390 ]; % 3000m
   %  % 
   %  % 
   %  depression_angle = zeros(size(depression_angles,1),30);
   %  for i = 1:size(depression_angles,1)
   %      depression_angle(i,:) = linspace(depression_angles(i,1),depression_angles(i,2),30);
   %  end

%% Calculate path loss between UAV and BS
%Path loss is depression angle dependant, Shadowing component is excluded
PL_UAV = zeros(size(d,2),size(depression_angle,2));
for i = 1:size(d,2)
    for j = 1:size(depression_angle,2)
        tmp1 = 10*alfa*log10(d(i)); % terrestrial path loss
        tmp2 = A*(depression_angle(i,j)-teta_0) * exp(-((depression_angle(i,j) - teta_0)/B)) + eta_0; % Aerial excess path-loss
       % tmp3 = normrnd(0,a*depression_angle(i,j)+sigma_0); Gaussian random variable - shadowing component.
        PL_UAV(i,j) = tmp1 + tmp2;% + tmp3;
    end
end  

%% Calculate SNR at UAV and convert it from dBm's to mW's
SNR_dBm = zeros(size(PL_UAV,1),size(PL_UAV,2));
SNR_mW = zeros(size(PL_UAV,1),size(PL_UAV,2));
Pt = 14; % dBm Transmit power
B = 180000; % Bandwidth
NF = 5; % Noise figure, https://altair.sony-semicon.com/wp-content/uploads/2017/02/Coverage-Analysis-of-LTE-CAT-M1-White-Paper.pdf
Pn = -174 + 10*log10(B) + NF; % Effective noise Figure
G_tx = 0;
G_rx = 0;


 for i = 1 : size(PL_UAV,1)
    for j = 1 : size(PL_UAV,2)
        SNR_dBm(i,j) = Pt - Pn - PL_UAV(i,j)  + G_rx + G_tx; % SNR in dBm
        SNR_mW(i,j) = 10^(SNR_dBm(i,j)/10); % Conversion to mW
    end
 end

B;
%% Helpers

%for i = 1 : length(d)
%    for j = 1 : size(depression_angle,2)
%        c = d(i)/cosd(depression_angle(i,j));
%        uav_altitude(i,j) = sqrt(c^2 - d(i)^2);
%    end
%end   
%     min = find(uav_altitude(i,:) > 15 & uav_altitude(i,:) < 16);
%     min = depression_angle(min(1));
%     max = find(uav_altitude(i,:) > 29 & uav_altitude(i,:) < 30)
%     if (depression_angle(max(1)) > 0 )
%         max = depression_angle(max(1))
%     else
%         max = -depression_angle(max(1))
%     end
%     uav_altitude_angles(i,:) = linspace(min,max,length(depression_angle));    
% end
% 
% for i = 1 : length(d)
%     for j = 1 : length(uav_altitude_angles)
%         c = d(i)/cosd(uav_altitude_angles(i,j));
%         uav_altitude_(i,j) = sqrt(c^2 - d(i)^2);
%     end
% end
