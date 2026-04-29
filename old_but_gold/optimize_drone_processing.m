function [local_opt, edge_opt, max_success] = optimal_packet_allocation(M, h_uav, SNR_dB, T, E, B)
%OPTIMAL_PACKET_ALLOCATION Determines the optimal allocation of M packets
% between local processing and edge transmission for a drone, maximizing
% successful processing under a battery constraint.
%
% Inputs:
%   M       - Total number of decoded packets (integer)
%   h_uav   - Drone altitude in meters (scalar)
%   SNR_dB  - Signal-to-noise ratio in dB (scalar)
%   T       - Predicted throughput (packets per slot, scalar)
%   E       - Predicted error rate on backhaul link (0 <= E < 1, scalar)
%   B       - Available battery energy in Joules (scalar)
%
% Outputs:
%   local_opt   - Optimal number of packets to process locally
%   edge_opt    - Optimal number of packets to transmit to edge server
%   max_success - Maximum number of successfully processed packets
%
% Example usage:
%   [local_i, edge_i, success] = optimal_packet_allocation(10, 100, 20, 8, 0.1, 5);

    %-------------------------------
    % Define system constants
    %-------------------------------
    E_local = 0.5;           % Energy per local processing (Joules/packet)
    PACKET_SIZE = 1000 * 8;  % Packet size in bits (1000 bytes)
    BANDWIDTH = 180e3;         % Channel bandwidth (Hz)
    FREQ = 2.4e9;            % Carrier frequency (Hz)
    N0 = 1e-20;              % Noise power spectral density (W/Hz)
    C = 3e8;                 % Speed of light (m/s)
    G_tx = 1;                % Transmit antenna gain (unitless)
    G_rx = 1;                % Receive antenna gain (unitless)

    % Convert SNR from dB to linear
    SNR_linear = 10^(SNR_dB / 10);

    % Precompute path loss (FSPL)
    d = h_uav; % Assume vertical link for simplicity
    lambda = C / FREQ;
    FSPL = (4 * pi * d / lambda)^2; % Linear scale

    % Precompute achievable data rate (Shannon capacity)
    R = BANDWIDTH * log2(1 + SNR_linear); % bits/sec

    % Precompute transmission time per packet
    t_tx = PACKET_SIZE / R; % seconds

    % Precompute required receive power for given SNR
    P_rx = SNR_linear * N0 * BANDWIDTH; % Watts

    % Precompute required transmit power
    P_tx = P_rx * FSPL / (G_tx * G_rx); % Watts

    % Transmission energy per packet
    E_tx = P_tx * t_tx; % Joules

    %-------------------------------
    % Initialize search
    %-------------------------------
    max_success = -Inf;
    local_opt = 0;
    edge_opt = 0;

    % Throughput constraint: maximum packets that can be transmitted
    max_edge = min(M, T);

    % Iterate over all feasible allocations
    for local_i = 0:M
        edge_i = M - local_i;

        % Enforce throughput constraint
        if edge_i > max_edge
            continue;
        end

        % Compute total energy consumption
        total_energy = local_i * E_local + edge_i * E_tx;

        % Check battery constraint
        if total_energy > B
            continue;
        end

        % Compute number of successfully processed packets
        success = local_i + edge_i * (1 - E);

        % Update optimal allocation if better
        if success > max_success
            max_success = success;
            local_opt = local_i;
            edge_opt = edge_i;
        end
    end

    % If no feasible allocation found, return zeros
    if max_success == -Inf
        local_opt = 0;
        edge_opt = 0;
        max_success = 0;
    end

    % Display results (optional)
    fprintf('Optimal allocation: %d local, %d edge\n', local_opt, edge_opt);
    fprintf('Maximum successful packets: %.2f\n', max_success);
    fprintf('Total energy used: %.3f J (Battery: %.3f J)\n', ...
        local_opt * E_local + edge_opt * E_tx, B);

end