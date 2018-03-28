function VAR = variance_ensemble(sys, N_it, T, dt, show_percentage)

% Inputs: state space, number of realizations, max time of one realization,
%         dt: sample rate of one realization

x = zeros(1,5);
var_data = zeros(N_it, 5);
rng('shuffle')
seed_i = randi([1 2^31],1,N_it);

for i = 1:N_it
    [~, V, a, t, q, N] = time_simulation(sys, dt, T, seed_i(i));
    var_data(i, :) = var([V, a, t, q, N], 1);
    perc = 100*i/N_it;
    if mod(perc, 1)==0 && show_percentage
        disp('% Done:')
        disp(perc)
    end
end

for i = 1:5
   x(i) = mean(var_data(:,i)); 
end

VAR = table(x(:, 1), x(:, 2), x(:, 3), x(:, 4), x(:, 5),...
    'VariableNames', {'sigma2_V','sigma2_alpha','sigma2_theta','sigma2_q','sigma2_Nz'});

end