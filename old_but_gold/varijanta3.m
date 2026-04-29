clear; clc;

%% 1. Inputs & System Parameters
M = 100; B = 30; T_max = 400;
E_cpu = 0.5; E_tx = 0.2;
L_loc = 10; L_edge = 5;

P_e_vec = [0.0086; 0.0682; 0.0368; 0.0271; 0.0307; 0.0466; 0.0790; 0.1377; ...
           0.2337; 0.3700; 0.5346; 0.7031; 0.8423; 0.9330; 0.9776; 0.9945];
num_scenarios = length(P_e_vec);

opt_p = zeros(num_scenarios, 1);
opt_q = zeros(num_scenarios, 1);

%% 2. Optimization Loop
for i = 1:num_scenarios
    Pe = P_e_vec(i);
    
    % f = Objective coefficients (We minimize -Success to maximize Success)
    % Variable vector x = [p; q]
    f = [-1; -(1 - Pe)];
    
    % A*x <= b (Inequality constraints)
    % Row 1: p + q <= M
    % Row 2: p*E_cpu + q*E_tx <= B
    % Row 3: p*L_loc <= T_max
    % Row 4: q*L_edge <= T_max
    A = [1, 1; 
         E_cpu, E_tx;
         L_loc, 0;
         0, L_edge];
    
    b = [M; B; T_max; T_max];
    
    % Integer constraints (both variables 1 and 2 are integers)
    intcon = [1, 2];
    
    % Lower bounds (p,q >= 0)
    lb = [0; 0];
    
    % Solve using intlinprog
    options = optimoptions('intlinprog','Display','off');
    [x, fval] = intlinprog(f, intcon, A, b, [], [], lb, [], [], options);
    
    opt_p(i) = x(1);
    opt_q(i) = x(2);
end

%% 3. Final Calculations & Display
success = opt_p + opt_q .* (1 - P_e_vec);
batt_rem = B - (opt_p * E_cpu + opt_q * E_tx);
time_used = T_max - max(opt_p * L_loc, opt_q * L_edge);

Results = table(P_e_vec, opt_p, opt_q, success, batt_rem, time_used, ...
    'VariableNames', {'Pe', 'Local_p', 'Edge_q', 'Success', 'BattRem', 'Time_left_ms'});
disp(Results);