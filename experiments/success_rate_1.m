clear all
close all
clc

addpath('./tools/')

%% test for network size, training length, 
% reset time and noise level
% please record the rmse and running time

% network_set = round(exp(linspace(log(20), log(270), 10)));
% training_length_set = round(exp(linspace(log(2000), log(150000), 10)));

reset_t_set = round(linspace(10, 150, 10));
noise_level_set = 10 .^ linspace(-3, 0, 10);

time_today = datestr(now, 'mmddyyyy');


%% first to test the n_set and train_len_set

% reset_t = 80;
% noise_level = 2.0 * 10 ^ (-2);
n = 100;
train_t = 10000;

bias = 2.0;

% delete(gcp('nocreate'))
% parpool('local',6)

iteration = 50;

rmse_set_lorenz = zeros(length(reset_t_set), length(noise_level_set), iteration);
rmse_set_circle = zeros(length(reset_t_set), length(noise_level_set), iteration);
rmse_set_mg17 = zeros(length(reset_t_set), length(noise_level_set), iteration);
rmse_set_infty = zeros(length(reset_t_set), length(noise_level_set), iteration);
time_set = zeros(length(reset_t_set), length(noise_level_set), iteration);

for rts = 1:length(reset_t_set)
    reset_t = reset_t_set(rts);
    for nls = 1:length(noise_level_set)
        noise_level = noise_level_set(nls);
        
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
        rmse_set_lorenz(rts, nls, :) = rmse_parfor_set_lorenz;
        rmse_set_circle(rts, nls, :) = rmse_parfor_set_circle;
        rmse_set_mg17(rts, nls, :) = rmse_parfor_set_mg17;
        rmse_set_infty(rts, nls, :) = rmse_parfor_set_infty;
        
        time_set(rts, nls, :) = time_parfor_set;
    end
end

save_success_rate.reset_t_set = reset_t_set;
save_success_rate.noise_level_set = noise_level_set;
save_success_rate.n = n;
save_success_rate.train_t = train_t;
save_success_rate.rmse_set_lorenz = rmse_set_lorenz;
save_success_rate.rmse_set_circle = rmse_set_circle;
save_success_rate.rmse_set_mg17 = rmse_set_mg17;
save_success_rate.rmse_set_infty = rmse_set_infty;
save_success_rate.time_set = time_set;

save(['save_data/save_success_rate_rn_', time_today, '_' num2str(randi(999)) '.mat'], 'save_success_rate')





