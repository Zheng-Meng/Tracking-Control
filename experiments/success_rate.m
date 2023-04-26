clear all
close all
clc

addpath('./tools/')

%% test for network size, training length, 
% reset time and noise level
% please record the rmse and running time

% network_set = round(linspace(20, 500, 10));
% training_length_set = round(linspace(20000, 200000, 10));

network_set = round(exp(linspace(log(20), log(250), 10)));
training_length_set = round(exp(linspace(log(2000), log(150000), 10)));

% network_set = fliplr(network_set);
% training_length_set = fliplr(training_length_set);

reset_t_set = round(linspace(10, 150, 10));
noise_level_set = 10 .^ linspace(-3, 0, 10);

time_today = datestr(now, 'mmddyyyy');


%% first to test the n_set and train_len_set

reset_t = 80;
noise_level = 2.0 * 10 ^ (-2);

bias = 2.0;

% delete(gcp('nocreate'))
% parpool('local',10)

iteration = 50;

rmse_set_lorenz = zeros(length(network_set), length(training_length_set), iteration);
rmse_set_circle = zeros(length(network_set), length(training_length_set), iteration);
rmse_set_mg17 = zeros(length(network_set), length(training_length_set), iteration);
rmse_set_infty = zeros(length(network_set), length(training_length_set), iteration);
time_set = zeros(length(network_set), length(training_length_set), iteration);

for ns = 1:length(network_set)
    n = network_set(ns);
    for tls = 1:length(training_length_set)
        train_t = training_length_set(tls);
        
        rmse_parfor_set_lorenz = zeros(1, iteration);
        rmse_parfor_set_circle = zeros(1, iteration);
        rmse_parfor_set_mg17 = zeros(1, iteration);
        rmse_parfor_set_infty = zeros(1, iteration);
        
        time_parfor_set = zeros(1, iteration);
        for repeat_i = 1:iteration
            [rmse_l, rmse_c, rmse_m, rmse_i, t_repeat_i] = func_train_val(n, train_t, reset_t, noise_level, bias);
            rmse_parfor_set_lorenz(repeat_i) = rmse_l;
            rmse_parfor_set_circle(repeat_i) = rmse_c;
            rmse_parfor_set_mg17(repeat_i) = rmse_m;
            rmse_parfor_set_infty(repeat_i) = rmse_i;
            
            time_parfor_set(repeat_i) = t_repeat_i;
        end

        rmse_set_lorenz(ns, tls, :) = rmse_parfor_set_lorenz;
        rmse_set_circle(ns, tls, :) = rmse_parfor_set_circle;
        rmse_set_mg17(ns, tls, :) = rmse_parfor_set_mg17;
        rmse_set_infty(ns, tls, :) = rmse_parfor_set_infty;
        
        time_set(ns, tls, :) = time_parfor_set;
    end
end

save_success_rate.reset_t = reset_t;
save_success_rate.noise_level = noise_level;
save_success_rate.network_set = network_set;
save_success_rate.training_length_set = training_length_set;
save_success_rate.rmse_set_lorenz = rmse_set_lorenz;
save_success_rate.rmse_set_circle = rmse_set_circle;
save_success_rate.rmse_set_mg17 = rmse_set_mg17;
save_success_rate.rmse_set_infty = rmse_set_infty;
save_success_rate.time_set = time_set;

save(['save_data/save_success_rate_nt_', time_today, '_' num2str(randi(999)) '.mat'], 'save_success_rate')




















































