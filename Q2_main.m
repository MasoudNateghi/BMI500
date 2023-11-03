clear; close all; clc; 
% parameters
beta = 0.3; 
gamma = 0.1;
sigma = 0.2; 
mu = 0.01; 

h = 1e-3; % step size
d = 365; % number of days
days = 1:h:d;

% initial values
S_0 = 990; 
E_0 = 9; 
I_0 = 1; 
R_0 = 0;

[S, E, I, R] = SEIR(beta, gamma, sigma, mu, d, h, S_0, E_0, I_0, R_0);

figure;
plot(days, S, 'LineWidth', 2)
hold on
plot(days, I, 'LineWidth', 2)
plot(days, R, 'LineWidth', 2)
plot(days, E, 'LineWidth', 2)
xlabel('days')
ylabel('Population Ratio')
legend('S(t)', 'I(t)', 'R(t)', 'E(t)', 'Interpreter', 'latex')
