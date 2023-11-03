clear; close all; clc; 
% parameters
beta = 0.3; 
gamma = 0.1;

h = 1e-5; % step size
d = 150; % number of days
days = 1:h:d; % create days (time) vector

% initial values
S_0 = 999;
I_0 = 1; 
R_0 = 0; 

[S, I, R] = SIR(beta, gamma, d, h, S_0, I_0, R_0);

figure;
plot(days, S, 'LineWidth', 2)
hold on
plot(days, I, 'LineWidth', 2)
plot(days, R, 'LineWidth', 2)
xlabel('days')
ylabel('Population Ratio')
legend('S(t)', 'I(t)', 'R(t)', 'Interpreter', 'latex')