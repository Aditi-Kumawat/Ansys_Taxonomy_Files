% Single Degree of Freedom System Response under Ricker Wavelet
clear; clc;

% System parameters
m = 1;       % Mass [kg]
k = 10;      % Stiffness [N/m]
c = 0.5;     % Damping coefficient [Ns/m]

% Ricker wavelet parameters
tMax = 5;    % Maximum time [s]
dt = 0.01;   % Time step [s]
f0 = 5;      % Peak frequency [Hz]

% Monte Carlo simulation parameters
numSimulations = 100;  % Number of simulations

% Generate time vector
t = 0:dt:tMax;

% Initialize response matrix
responseMatrix = zeros(length(t), numSimulations);

% Perform Monte Carlo simulations
for i = 1:numSimulations
    % Generate random initial conditions
    x0 = randn(1);      % Initial displacement
    v0 = randn(1);      % Initial velocity
    
    % Calculate natural frequency and critical damping coefficient
    wn = sqrt(k/m);
    cc = 2 * sqrt(k*m);
    
    % Initialize response vectors
    x = zeros(size(t));
    v = zeros(size(t));
    
    % Apply Ricker wavelet
    wavelet = (1 - 2 * pi^2 * f0^2 * (t - tMax/2).^2) .* exp(-pi^2 * f0^2 * (t - tMax/2).^2);
    
    % Perform time integration using Newmark's method
    for j = 2:length(t)
        dt = t(j) - t(j-1);
        
        % Calculate acceleration using the equation of motion
        a = (1/m) * (wavelet(j) - c*v(j-1) - k*x(j-1));
        
        % Update displacement and velocity
        x(j) = x(j-1) + dt*v(j-1) + (dt^2/2)*(1-2*0.6)*a;
        v(j) = v(j-1) + dt*((1-0.6)*a + 0.6*((1-0.6)*a + (1/m)*wavelet(j)));
    end
    
    % Store the response in the response matrix
    responseMatrix(:, i) = x';
end

% Plot the Monte Carlo responses
figure;
for i = 1:numSimulations
    plot(t, responseMatrix(:, i));
    hold on;
end
xlabel('Time [s]');
ylabel('Displacement [m]');
title('Monte Carlo Simulation: Single Degree of Freedom System Response');
grid on;
hold off;
