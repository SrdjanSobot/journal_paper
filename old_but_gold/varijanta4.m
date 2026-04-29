clear; clc;

%% 1. Parameters
B_current = 20;     % Starting Battery (Joules)
E_cpu = 0.5; E_tx = 0.2;
L_loc = 10; L_edge = 5;
T_max = 400;        % Slot deadline
lambda = 1.7; 

% (m packets decoded per slot, with cumsum 100 )
rng(1); % For reproducibility
M_arrivals = randsample(1:3, 100, true);

% Original Vector
P_e_vec = [0.0086, 0.0682, 0.0368, 0.0271, 0.0307, 0.0466, 0.0790, 0.1377, ...
           0.2337, 0.3700, 0.5346, 0.7031, 0.8423, 0.9330, 0.9776, 0.9945];

P_e_vec_100 = interp1(1:16, P_e_vec, linspace(1, 16, 100), 'linear');


results = zeros(length(M_arrivals), 6);

%% 2. Real-Time Processing Loop
for t = 1:length(M_arrivals)
    if B_current <= 0, break; end
    
    num_packets_in_slot_t = M_arrivals(t);
    Pe = P_e_vec_100(t);
    
    % Calculate Current Efficiencies
    %eff_loc = 1.0 / E_cpu;            % Always 2.0 for your parameters
    %eff_edge = (1 - Pe) / E_tx;       % Variable based on Pe
    U_loc  = 1.0 - (lambda * E_cpu);
    U_edge = (1 - Pe) - (lambda * E_tx);
    % Optimization: Pick the best mode that beats the threshold
    if U_edge >= U_loc && U_edge > 0
        % Edge is better
        edge = min(num_packets_in_slot_t, floor(T_max/L_edge));
        loc = 0;
    elseif U_loc > U_edge && U_loc > 0
        % Local is better and worth the energy
        loc = min(num_packets_in_slot_t, floor(T_max/L_loc));
        edge = 0;
    else
        % Channel is too bad/Energy too expensive: Drop packets to save mission
        loc = 0; edge = 0;
    end
    
    % Update Battery and Success
    cost = (loc * E_cpu) + (edge * E_tx);
    B_current = B_current - cost;
    success = loc + edge*(1-Pe);
    
    results(t, :) = [Pe, num_packets_in_slot_t, loc, edge, success, B_current];
end

%% 3. Display
% T = table(results(:,1), results(:,2), results(:,3), results(:,4), results(:,5), results(:,6), ...
%     'VariableNames', {'Pe', 'Arrived', 'Local_p', 'Edge_q', 'Success', 'BattRem'});
% disp(T);

%% 3. Summary Calculations
total_arrived = sum(results(:, 2));
total_local   = sum(results(:, 3));
total_edge    = sum(results(:, 4));
total_success = sum(results(:, 5));
total_dropped = total_arrived - (total_local + total_edge);
final_battery = B_current;
success_rate  = (total_success / total_arrived) * 100;

%% 4. Display Summary Report
fprintf('==========================================\n');
fprintf('            MISSION SUMMARY               \n');
fprintf('==========================================\n');
fprintf('Total Packets Collected (Arrived): %d\n', total_arrived);
fprintf('Total Packets Processed Locally:   %d\n', total_local);
fprintf('Total Packets Offloaded to Edge:   %d\n', total_edge);
fprintf('Total Packets Dropped (Battery):   %d\n', total_dropped);
fprintf('------------------------------------------\n');
fprintf('TOTAL SUCCESSFUL PACKETS:          %.2f\n', total_success);
fprintf('OVERALL SUCCESS RATE:              %.2f%%\n', success_rate);
fprintf('REMAINING BATTERY:                 %.2f J\n', final_battery);
fprintf('==========================================\n');