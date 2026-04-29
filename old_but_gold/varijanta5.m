clear; clc;

%% 1. Parameters & Setup
B_start = 20; 
E_cpu = 0.5; E_tx = 0.2;
L_loc = 10; L_edge = 5; T_max = 400;

% Adaptive Lambda Range
lambda_min = 0.8;   % Aggressive (Full battery)
lambda_max = 4.0;   % Conservative (Low battery)
k_sensitivity = 2;  % Quadratic curve

% Up-sampled P_e (100 points)
P_e_orig = [0.0086, 0.0682, 0.0368, 0.0271, 0.0307, 0.0466, 0.0790, 0.1377, ...
            0.2337, 0.3700, 0.5346, 0.7031, 0.8423, 0.9330, 0.9776, 0.9945];
P_e_vec_100 = interp1(1:16, P_e_orig, linspace(1, 16, 100), 'linear');
rng(1); % For reproducibility
M_arrivals = randsample(1:3, 100, true); 

%% 2. The Decision Loop
B_curr = B_start;
results = zeros(100, 7); % Added a column for lambda

for t = 1:100
    if B_curr <= 0, break; end
    
    % --- ADAPTIVE STEP ---
    % Calculate lambda based on current battery percentage
    batt_ratio = B_curr / B_start;
    current_lambda = lambda_min + (lambda_max - lambda_min) * (1 - batt_ratio)^k_sensitivity;
    
    m = M_arrivals(t);
    Pe = P_e_vec_100(t);
    
    % Utility Calculation
    U_loc  = 1.0 - (current_lambda * E_cpu);
    U_edge = (1 - Pe) - (current_lambda * E_tx);
    
    % Decision Logic
    if U_edge >= U_loc && U_edge > 0
        q = min(m, floor(T_max/L_edge)); p = 0;
    elseif U_loc > U_edge && U_loc > 0
        p = min(m, floor(T_max/L_loc)); q = 0;
    else
        p = 0; q = 0; % Drop
    end
    
    % Update Stats
    cost = (p * E_cpu) + (q * E_tx);
    B_curr = B_curr - cost;
    success = p + q*(1 - Pe);
    
    % Log: [Pe, m, p, q, success, B_curr, current_lambda]
    results(t, :) = [Pe, m, p, q, success, B_curr, current_lambda];
end

%% 3. Visualization of Adaptivity
figure('Color', 'w');
yyaxis left
plot(results(:, 6), 'LineWidth', 2); ylabel('Remaining Battery (J)');
hold on;
yyaxis right
plot(results(:, 7), 'r--', 'LineWidth', 2); ylabel('Adaptive \lambda (Penalty)');
xlabel('Time Slot'); grid on;
title('Adaptive Energy Awareness vs. Battery Depletion');
legend('Battery Level', 'Adaptive \lambda Value');

%% 4. Summary Calculation
total_arrived = sum(results(:, 2));
total_success = sum(results(:, 5));
fprintf('Total Success: %.2f packets\n', total_success);
fprintf('Success Rate: %.2f%%\n', (total_success/total_arrived)*100);

%% 3. Summary Calculations
total_arrived = sum(results(:, 2));
total_local   = sum(results(:, 3));
total_edge    = sum(results(:, 4));
total_success = sum(results(:, 5));
total_dropped = total_arrived - (total_local + total_edge);
final_battery = B_curr;
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