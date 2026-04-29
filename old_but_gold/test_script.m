% Iterate through all possible local processing amounts (0 to M)
clear;clc;
M = 50;
E = load("Pe_test.mat","P_e");
E = E.P_e(:,1);
max_processed = zeros(length(E),M)*-inf;
best_energy = zeros(length(E),M)*-inf;
valid_solution_found(length(E),M) = false;   
local = zeros(length(E),M)*-inf;
edge = zeros(length(E),M)*-inf;
E_cpu_per_packet = 0.5;
E_tx_per_packet = 0.2;

B = 15;
for i = 1:length(E)
    current_energy = 0;
    current_success = 0;
    for p = 0:M
        q = M - p; % The rest go to the edge
        
        % Calculate total energy required for this split
        current_energy = (p * E_cpu_per_packet) + (q * E_tx_per_packet);
        
        % Check if this split is feasible within Battery limits
        if current_energy <= B
            valid_solution_found(i,p+1) = true;
            
            % Calculate Expected Successful Packets
            % Local packets have 100% success (1.0). Edge packets have (1-E) success.
            current_success = (p * 1.0) + (q * (1 - E(i)));
            
            max_processed(i, p+1) = current_success;
            local(i, p+1) = p;
            edge(i, p+1) = q;
            best_energy(i, p+1) = current_energy;
        end
    end
end