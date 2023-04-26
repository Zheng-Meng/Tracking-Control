clear all
close all
clc

addpath('./tools/')

% load training data
load('./save_file/all_traj_06282022.mat')
time_today = datestr(now, 'mmddyyyy');

plot_movie = 0;
if exist('traj_type','var') == 0
    traj_type = 'infty';
end

traj_set = ["lorenz", "circle", "mg17", "infty", "astroid", "fermat", ...
    "lissajous", "talbot", "heart", "chua", "rossler", ...
    "sprott_1", "sprott_4", "mg30", "epitrochoid"];
% traj_set = ["lorenz"];
iteration = 100;

bridge_type = 'cubic';
time_infor.val_length=250000;
val_length_rm = time_infor.val_length;

failure.type = 'none';
blur.last_time = 1000;
blur.recover_time = time_infor.val_length;
blur.blur = 0;

failure.type = 'all';
failure.amplitude = 0.1;
failure.amplitude_2 = 0.1;

rmse_start_time = round(val_length_rm * 3/5);
rmse_end_time = time_infor.val_length - 100;

rmse_set = zeros(length(traj_set), iteration);

idx=1;
for traj_id = 1:length(traj_set)
    traj_type = traj_set(traj_id);
    if strcmp(traj_type, 'circle') == 1
        traj_frequency = 150;
    elseif strcmp(traj_type, 'infty') == 1
        traj_frequency = 75;
    end

    rmse_parfor_set = zeros(1, iteration);
    for repeat_i = 1:iteration
        load('./save_file/all_traj_06282022.mat')
        time_infor.val_length = val_length_rm;
        save_rend=0;
        val_and_update;
        rmse_parfor_set(repeat_i) = func_rmse(data_pred, data_control, rmse_start_time, rmse_end_time);
        aaa = 1;
    end
    rmse_parfor_set = sort(rmse_parfor_set);
    rmse_set(traj_id, :) = rmse_parfor_set;
end

save_success.rmse_set = rmse_set;
save_success.val_length = time_infor.val_length;
save_success.traj_set = traj_set;

% save(['./save_data/save_saferegion_success_rate_noise_' time_today, '_' num2str(randi(999)) '.mat'], "save_success")

%%
% load('./save_data/save_normal_success_rate_08022022_92.mat')

% rmse_set = save_success.rmse_set;

rmse_threshold = 0.18;
rmse_logic = zeros(size(rmse_set));
for i = 1:size(rmse_set, 1)
    for j = 1:size(rmse_set, 2)
        if rmse_set(i, j) > rmse_threshold
            rmse_logic(i, j) = 0;
        else
            rmse_logic(i, j) = 1;
        end
    end
end

rmse_count = mean(rmse_logic, 2);

figure()
plot(rmse_count, 'o')




















