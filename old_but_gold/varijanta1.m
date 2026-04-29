clear; clc;

% --- Inputs ---
M = 100;
B = 30; 
E_cpu = 0.5;
E_tx = 0.2;

% Simulated Pe data based on your input
P_e_vec = [0.0086; 0.0682; 0.0368; 0.0271; 0.0307; 0.0466; 0.0790; 0.1377; ...
           0.2337; 0.3700; 0.5346; 0.7031; 0.8423; 0.9330; 0.9776; 0.9945];

num_scenarios = length(P_e_vec);

% --- 1. Calculate Efficiencies ---
eff_local = 1.0 / E_cpu;            % Success per Joule (Local)
eff_edge = (1 - P_e_vec) / E_tx;    % Success per Joule (Edge) - Vector

% --- 2. Optimization Logic (Vectorized) ---
% We create placeholders for our decisions
p = zeros(num_scenarios, 1);
q = zeros(num_scenarios, 1);

% Identify where Edge is more efficient than Local
edge_is_better = eff_edge >= eff_local;

% CASE A: Edge is more efficient (Prioritize q)
idx_e = edge_is_better;
if any(idx_e)
    % Maximize q based on battery
    max_q_budget = floor(B / E_tx);
    q(idx_e) = min(M, max_q_budget);
    
    % If battery remains, "Upgrade" edge packets to local for 100% success
    energy_left = B - (q(idx_e) * E_tx);
    swap_cost = E_cpu - E_tx;
    possible_swaps = floor(energy_left / swap_cost);
    
    % Number of packets to move from Edge to Local
    num_to_upgrade = min(q(idx_e), possible_swaps);
    
    p(idx_e) = num_to_upgrade;
    q(idx_e) = q(idx_e) - num_to_upgrade;
end

% CASE B: Local is more efficient (Prioritize p)
idx_l = ~edge_is_better;
if any(idx_l)
    % Maximize p based on battery
    max_p_budget = floor(B / E_cpu);
    p(idx_l) = min(M, max_p_budget);
    
    % If battery remains, use it for Edge processing
    energy_left = B - (p(idx_l) * E_cpu);
    q(idx_l) = min(M - p(idx_l), floor(energy_left / E_tx));
end

% --- 3. Results ---
expected_success = (p * 1.0) + (q .* (1 - P_e_vec));

T = table(P_e_vec, p, q, expected_success, ...
    'VariableNames', {'ErrorRate', 'Local_p', 'Edge_q', 'Success'});
disp(T);