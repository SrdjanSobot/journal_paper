clear;
% 1. System Constants
E_tx = 0.2; E_cpu = 0.5;
lambda = 1.5; % The "Mission Duration" penalty
P_e = load('Pe_test.mat', 'P_e');
P_e_matrix = P_e.P_e;

results = [];

for i = 1:size(P_e_matrix, 1)
    current_Pe = P_e_matrix(i, :);
    
    % Step 1: Calculate Utility for all 5 Edge Packing Options
    % Utility = (Expected Success) - (Penalty * Energy)
    % Expected Success = j * (1 - Pe)
    j_options = 1:5;
    U_edge_options = j_options .* (1 - current_Pe) - (lambda * E_tx);
    
    % Step 2: Calculate Local Utility (per single packet)
    U_local_single = 1.0 - (lambda * E_cpu);
    
    % Step 3: Find the "Winner"
    [max_U_edge, best_j] = max(U_edge_options);
    
    if max_U_edge > U_local_single && max_U_edge > 0
        decision = sprintf('Edge (Pack %d)', best_j);
        success = best_j * (1 - current_Pe(best_j));
        energy = E_tx;
    elseif U_local_single > 0
        decision = 'Local';
        success = 1;
        energy = E_cpu;
    else
        decision = 'Dropped';
        success = 0;
        energy = 0;
    end
    
    results = [results; i, best_j, success, energy];
    fprintf('Scenario %d: Decision = %s | Expected Success = %.4f\n', i, decision, success);
end