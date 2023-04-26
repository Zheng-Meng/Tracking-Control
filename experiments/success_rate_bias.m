clear all
close all
clc

addpath('./tools/')

bias_set = linspace(0, 3, 7);

time_today = datestr(now, 'mmddyyyy');

reset_t = 80;
noise_level = 2.0 * 10 ^ (-2);
n=200;
train_t = 150000;

iteration = 10;

rmse_set_lorenz = zeros(length(bias_set), iteration);
rmse_set_circle = zeros(length(bias_set), iteration);
rmse_set_mg17 = zeros(length(bias_set), iteration);
rmse_set_infty = zeros(length(bias_set), iteration);

for bs = 1:length(bias_set)
    bias = bias_set(bs);
    rmse_parfor_set_lorenz = zeros(1, iteration);
    rmse_parfor_set_circle = zeros(1, iteration);
    rmse_parfor_set_mg17 = zeros(1, iteration);
    rmse_parfor_set_infty = zeros(1, iteration);
    
    for repeat_i = 1:iteration
        [rmse_l, rmse_c, rmse_m, rmse_i, t_repeat_i] = func_train_val(n, train_t, reset_t, noise_level, bias);
        rmse_parfor_set_lorenz(repeat_i) = rmse_l;
        rmse_parfor_set_circle(repeat_i) = rmse_c;
        rmse_parfor_set_mg17(repeat_i) = rmse_m;
        rmse_parfor_set_infty(repeat_i) = rmse_i;
    end
    
    rmse_set_lorenz(bs, :) = rmse_parfor_set_lorenz;
    rmse_set_circle(bs, :) = rmse_parfor_set_circle;
    rmse_set_mg17(bs, :) = rmse_parfor_set_mg17;
    rmse_set_infty(bs, :) = rmse_parfor_set_infty;
end

save_success_rate.reset_t = reset_t;
save_success_rate.noise_level = noise_level;
save_success_rate.n = n;
save_success_rate.train_t = train_t;
save_success_rate.bias_set = bias_set;
save_success_rate.rmse_set_lorenz = rmse_set_lorenz;
save_success_rate.rmse_set_circle = rmse_set_circle;
save_success_rate.rmse_set_mg17 = rmse_set_mg17;
save_success_rate.rmse_set_infty = rmse_set_infty;


save(['save_data/save_success_rate_bias_', time_today, '_' num2str(randi(999)) '.mat'], 'save_success_rate')





























