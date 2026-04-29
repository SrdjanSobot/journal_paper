clear; clc;

% --- Inputs ---
M = 100;
B = 1;
E_cpu_per_packet = 0.5;
E_tx_per_packet = 0.2;

% Load your error rates (simulating data for this example)
% Assuming E is a column vector of error rates
E_data = load("Pe_test.mat", "P_e");
P_e_vec = E_data.P_e(:,1); 
num_scenarios = length(P_e_vec);

% --- Optimization (Vectorized) ---

% 1. Determine max possible local packets given battery constraint
% Derivation: p*E_cpu + (M-p)*E_tx <= B
%             p*(E_cpu - E_tx) <= B - M*E_tx
if E_cpu_per_packet > E_tx_per_packet
    max_p_energy = floor((B - M * E_tx_per_packet) / (E_cpu_per_packet - E_tx_per_packet));
else
    % If local is cheaper than tx, we process EVERYTHING locally (p=M)
    max_p_energy = M; 
end

% 2. Enforce limits: p must be between 0 and M
optimal_p = min(M, max(0, max_p_energy));
optimal_q = M - optimal_p;

% 3. Check Feasibility
% Calculate energy used for this optimal split
energy_used = (optimal_p * E_cpu_per_packet) + (optimal_q * E_tx_per_packet);
is_feasible = energy_used <= B;

% 4. Calculate Objective (Expected Processed Packets)
% We replicate optimal_p to match the size of P_e_vec
optimal_p_vec = repmat(optimal_p, num_scenarios, 1);
optimal_q_vec = repmat(optimal_q, num_scenarios, 1);

% Success = Local + Edge_Success
% Note: Local is 100% success. Edge is (1 - P_e).
% Since p and q are constant regardless of P_e (because constraints 
% rely on Energy, not Error), p is the same for all error rates.
max_processed = (optimal_p_vec * 1.0) + (optimal_q_vec .* (1 - P_e_vec));

% If the configuration violates battery (negative energy budget), set result to NaN or 0
max_processed(~is_feasible) = NaN;

% --- Display Results for first 5 scenarios ---
T = table(P_e_vec(1:5), optimal_p_vec(1:5), optimal_q_vec(1:5), max_processed(1:5), ...
    'VariableNames', {'ErrorRate', 'Local_p', 'Edge_q', 'ExpectedSuccess'});
disp(T);