clear;
%T = load('T.mat', 'T');
%T = T.T;
%P_e = load('P_e.mat', 'P_e');
%P_e = P_e.P_e;
%Cmr = load('cmr.mat','Cmr');
%Cmr = Cmr.Cmr;

%% ######### Inputs & System Parameters #########
% Tier 1
M = 5; % MPR capability
MM = 11; % collect them all
U = [150 140 130 120 110 100 90 80 70 60]; % num of users
h_uav = [-20 40 60 80 100 120 140 160 180 200]; % UAV altitude
p_u = 0.1; % activation probability
b = 0.1;

% Tier 2
T_f = 100 * 10^-3; % Frame duration
R = 32 * 10^3; % Digital rate in kbps
n = T_f * R; % num of channel uses
k = 640;
RR = [0.6 0.7 0.8 0.9 1];
D_f = 250 * 10^3;

% Local resources
B = 25;                 % Battery budget (Joules)
T_max = 400;            % Latency Deadline (ms)

% Energy Costs (Joules per packet)
E_cpu = 0.5;            
E_tx = 0.2;             

% Latency Costs (ms per packet)
L_loc = 10;             % Time to process 1 packet locally
L_edge = 5;              % Time to transmit + process 1 packet at edge

% Limits
max_loc_time = floor(T_max / L_loc);   % Max local packets before deadline
max_edge_time = floor(T_max / L_edge);  % Max edge packets before deadline
lambda = 1.5; % Mission duration penality



%% Tier 1 calculations 
tic
if ~exist('Cmr', 'var')
    Cmr = zeros(length(h_uav), MM+1);
    for h = 1 : length(h_uav)
        for m = 0 : MM
             for i = 0 : U(h)
                 Cmr(h,m+1) = Cmr(h,m+1) + CalculateCmr_mex( int16(i), int16(m), h_uav(h), b ) * binopdf(i,U(h),p_u);
             end
        end
    end
    save('cmr.mat','Cmr');
end
toc

%% Tier 2 calculations
if ~exist('SNR', 'var')
    [h_uav, SNR, SNR_dBm] = snr_from_excess_PL_novo(h_uav);
    SNR = SNR';
end
packet_distribution = repmat(struct('edge',0,'local',0), size(SNR,1), size(SNR,2));
if ~exist('T', 'var')
    T = zeros(size(SNR,1), size(SNR,2));
    eff = zeros(size(SNR,1), size(SNR,2));
    P_e = zeros(size(SNR,2),M);
    for s_row = 1 : size(SNR,1)
        P_e = zeros(M, size(SNR,2));
        for s_column = 1 : size(SNR,2)
            for m = 1 : M       
                P_e(m, s_column) = frame_error_probability(SNR(s_row,s_column), n, m, k, 0);
            end
            current_Pe = P_e(:, s_column);
            eff_local = 1.0 / E_cpu;
            eff_edge = (1 - current_Pe) / E_tx;
            packet_distribution(s_row, s_column).edge = sum(eff_edge > eff_local*2);
            packet_distribution(s_row, s_column).local = sum(eff_edge < eff_local*2);

                % OVDE SAM STAO
                % ALGORITAM NA OSNOVU p_e ODLUČI KOLIKO ŠALJE NA EDGE
                % KOLIKO PROCESIRA LOKALNO.
                % DODATI DALJI TOK ALGORITMA, I PREBACITI SVE U OVU PETLJU.
                % KOLIKO JE OTIŠLO LOKALNO, KOLIKO NA EDGE, UPDATE STATUSA
                % BATERIJA, ZATIM ISPIS REZULTATA I TJT. PRETHODNI KOD
                % ISPOD VEEEROVATNO NIJE POTREBAN I SAV PROCESSING CE BITI
                % U OVOJ PETLJI
        end
    end
end
edge_matrix = arrayfun(@(x) x.edge, packet_distribution);
    save('T.mat','T');
    save('eff.mat','eff');
    save('P_e.mat','P_e');


%% Prepare for the impact!!! 
num_scenarios = size(P_e);
optimal_loc = zeros(num_scenarios);
optimal_edge = zeros(num_scenarios);

% Effectiveness
eff_local = 1.0 / E_cpu;
eff_edge = (1 - P_e) / E_tx;

%% CASE A: Edge is more efficient (Low Error Rate)
idx_e = eff_edge >= eff_local;
if any(idx_e)
    max_edge_battery = floor(B / E_tx);
    edge_vals = min([M, max_edge_time, max_edge_battery]);
    
    % Use remaining battery for Local processing
    battery_left = B - (edge_vals * E_tx);
    max_loc_battery = floor(battery_left / E_cpu);
    loc_vals = min([M - edge_vals, max_loc_time, max_loc_battery]);

     if(loc_vals == 0)
        swap = min([max_loc_time, max_loc_battery, M]);
        edge_vals = edge_vals - swap;
        loc_vals =  loc_vals + swap;
        optimal_edge(idx_e) = edge_vals;
        optimal_loc(idx_e) = loc_vals;
     else
        optimal_edge(idx_e) = edge_vals;
        optimal_loc(idx_e) = loc_vals;
     end
    battery_left = B - (loc_vals * E_cpu + edge_vals * E_tx);
    local_time_left = T_max - (loc_vals * L_loc); % Parallel model
    
    while((battery_left >= E_cpu && local_time_left > 0))
        max_loc_battery = floor(battery_left / E_cpu);
        swap = min(local_time_left, max_loc_battery);
        if( swap > loc_vals || edge_vals < swap )
            break;
        end
        edge_vals = edge_vals - swap;
        loc_vals =  loc_vals + swap;
        optimal_edge(idx_e) = edge_vals;
        optimal_loc(idx_e) = loc_vals;
        battery_left = B - (loc_vals * E_cpu + edge_vals * E_tx);
        local_time_left = T_max - (loc_vals * L_loc); % Parallel model
        
    end
end

%% Case B: Local is more efficient (High error rate)
idx_l = ~idx_e;
if any(idx_l)
    max_loc_battery = floor(B / E_cpu);
    loc_vals = min([M, max_loc_time, max_loc_battery]);
    battery_left = B - (loc_vals * E_cpu);
    max_edge_battery = floor(battery_left / E_tx);
    edge_vals = min([M - loc_vals, max_edge_time, max_edge_battery]);
    
    optimal_loc(idx_l) = loc_vals;
    optimal_edge(idx_l) = edge_vals;
end

%% ---Final Calculations ---
success = (optimal_loc * 1.0) + (optimal_edge .* (1 - P_e));
batt_left = B - (optimal_loc * E_cpu + optimal_edge * E_tx);
time_used = max(optimal_loc * L_loc, optimal_edge * L_edge); % Parallel model

%% ---Display ---
T = table(P_e, optimal_loc, optimal_edge, success, batt_left, time_used, ...
    'VariableNames', {'Pe', 'Local_p', 'Edge_q', 'Success', 'BattRem', 'TimeMs'});
disp(T);


%% Helper functions
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

function [edge local] = split_or_not_to_split(P_e)
    eff_edge = (1 - P_e) / E_tx;
    mean(floor(eff_edge .* 100))
    eff_local = 1.0 / E_cpu;
    idx_e = eff_edge >= eff_local;
end

function batt_left = battery_status(loc, edge, B)
    batt_left = B - (loc * E_cpu + edge * E_tx);
end

