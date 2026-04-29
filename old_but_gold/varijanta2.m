clear;
T = load('C:\Users\srdja\Documents\MATLAB\simulacije\matlab\main_novo\slike_final\prva, snr_usecase1.mat', 'T');
T = T.T;
%% Inputs & System Parameters
M = 100;                % Total packets to process
B = 25;                 % Battery budget (Joules)
T_max = 400;            % Latency Deadline (ms)

% Energy Costs (Joules per packet)
E_cpu = 0.5;            
E_tx = 0.2;             

% Latency Costs (ms per packet)
L_loc = 10;             % Time to process 1 packet locally
L_edge = 5;              % Time to transmit + process 1 packet at edge

% error rates
P_e_vec = [0.0086; 0.0682; 0.0368; 0.0271; 0.0307; 0.0466; 0.0790; 0.1377; ...
           0.2337; 0.3700; 0.5346; 0.7031; 0.8423; 0.9330; 0.9776; 0.9945];

num_scenarios = length(P_e_vec);
optimal_loc = zeros(num_scenarios, 1);
optimal_edge = zeros(num_scenarios, 1);

%% Limits
max_loc_time = floor(T_max / L_loc);   % Max local packets before deadline
max_edge_time = floor(T_max / L_edge);  % Max edge packets before deadline

%
eff_local = 1.0 / E_cpu;
eff_edge = (1 - P_e_vec) / E_tx;

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
        swap = min(max_loc_time, max_loc_battery);
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
success = (optimal_loc * 1.0) + (optimal_edge .* (1 - P_e_vec));
batt_left = B - (optimal_loc * E_cpu + optimal_edge * E_tx);
time_used = max(optimal_loc * L_loc, optimal_edge * L_edge); % Parallel model

%% ---Display ---
T = table(P_e_vec, optimal_loc, optimal_edge, success, batt_left, time_used, ...
    'VariableNames', {'Pe', 'Local_p', 'Edge_q', 'Success', 'BattRem', 'TimeMs'});
disp(T);