function [S, I, R] = SIR(beta, gamma, d, h, S_0, I_0, R_0)
    days = 1:h:d;
    n_sample = length(days); % number of samples
    
    % create SIR vectors
    S = zeros(n_sample, 1); 
    I = zeros(n_sample, 1);
    R = zeros(n_sample, 1);
    
    % initial values
    S(1) = S_0; 
    I(1) = I_0; 
    R(1) = R_0; 
    
    
    % euler method to solve ODE
    for n = 1:length(days) - 1
        S(n + 1) = S(n) + h * (-beta * S(n) * I(n)); 
        I(n + 1) = I(n) + h * (beta * S(n) * I(n) - gamma * I(n));
        R(n + 1) = R(n) + h * gamma * I(n);
    end
end