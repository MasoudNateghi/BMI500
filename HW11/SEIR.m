function [S, E, I, R] = SEIR(beta, gamma, sigma, mu, d, h, S_0, E_0, I_0, R_0)
    days = 1:h:d;
    n_sample = length(days); % number of samples
    
    % create SIR vectors
    S = zeros(n_sample, 1); 
    E = zeros(n_sample, 1); 
    I = zeros(n_sample, 1);
    R = zeros(n_sample, 1);
    N = zeros(n_sample, 1); 
    
    % initial values
    S(1) = S_0; 
    E(1) = E_0;
    I(1) = I_0; 
    R(1) = R_0;
    
    % euler method to solve ODE
    for n = 1:length(days) - 1
        N(n) = S(n) + E(n) + I(n) + R(n);
        S(n + 1) = S(n) + h * (mu * N(n) - beta * S(n) * I(n) - mu * S(n)); 
        E(n + 1) = E(n) + h * (beta * S(n) * I(n) - (sigma + mu) * E(n));
        I(n + 1) = I(n) + h * (sigma * E(n) - (gamma + mu) * I(n));
        R(n + 1) = R(n) + h * (gamma * I(n) - mu * R(n));
    end
end