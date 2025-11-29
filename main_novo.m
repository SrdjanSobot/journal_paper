%% Tier 1
M = 3; % MPR capability
MM = 11; % collect them all
U = [150 140 130 120 110 100 90 80 70]; % num of users
L = [40 60 80 100 120 140 160 180 200]; % UAV altitude
p_u = 0.1; % activation probability
b = 0.1;

%% Tier 2
T_f = 100 * 10^-3; % Frame duration
R = 32 * 10^3; % Digital rate in kbps
n = T_f * R; % num of channel uses
k = 640;
RR = 0.6;
D_f = 250 * 10^3;

%% Tier 1 
if ~exist('SNR', 'var')
    [h_uav, SNR, SNR_dBm] = snr_from_excess_PL_novo(L);
    SNR = SNR';
end

T = zeros(size(SNR,1), size(SNR,2));
eff = zeros(size(SNR,1), size(SNR,2));
P_e = zeros(size(SNR,2),M);

tic
if ~exist('Cmr', 'var')
    Cmr = zeros(length(h_uav), MM+1);
    for h = 1:length(L)
        for m = 0 : MM
             for i = m : U(h)
                 Cmr(h,m+1) = Cmr(h,m+1) + CalculateCmr_mex( int16(i), int16(m), h_uav(h), b ) * binopdf(i,U,p_u);
             end
        end
    end
end
toc


%% Tier 2
for s_row = 1 : size(SNR,1)
    P_e = zeros(size(SNR,2),M);
    for s_column = 1 : size(SNR,2)
        for m = 1 : M

            P_e(s_column,m) = frame_error_probability(SNR(s_row,s_column), n, m, k, RR);
            
            if(m < M)
                T(s_row, s_column) = T(s_row, s_column) + ...
                                    (Cmr(s_row,m) * (m * (1 - P_e(s_column,m))));
            else
                T(s_row, s_column) = T(s_row, s_column) + ...
                                    (sum(Cmr(s_row,m:end)) * (m * (1 - P_e(s_column,m))));
            end
            if( m < M)
                eff(s_row, s_column) = eff(s_row, s_column) + efficiency(n, m, k, Cmr(s_row,m), RR);
            else
                eff(s_row, s_column) = eff(s_row, s_column) + efficiency(n, m, k, sum(Cmr(s_row,m:end)), RR);
            end
        end
    end
end

%% helper functions

function [P_e] = frame_error_probability(SNR, n, m, k, R)
    if R == 0
        R = k/n * m;
    end
    N_samples = 10^6;
    Nak_m = 3; % Nakagami paramter
    gRF = gamrnd(Nak_m,1/Nak_m,1,N_samples);
    SNR_Nak = SNR .* gRF;
    V_Nak = (1 - 1 ./ (1 + (SNR_Nak).^2)) .* (log2(exp(1))^2);
    C_Nak = log2(1+SNR_Nak);
    epsilon_Nak_all = qfunc(sqrt(n./V_Nak) .* (C_Nak-R));
    P_e = sum(epsilon_Nak_all) / N_samples;
end



function [eff] = efficiency (n, m, k, Cmr, R)
    if R == 0
        R = k/n * m;
    end
    
    eff = ( ( ( k * m ) / R ) * Cmr ) / n;
end
