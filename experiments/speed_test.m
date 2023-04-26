% clear all
% close all
% clc

addpath('./tools/')

%%

load('./save_file/all_traj_06282022.mat')

time_today = datestr(now, 'mmddyyyy');

plot_val_and_update = 0;

blur.blur = 0;
plot_movie = 0;
traj_type = 'infty';

if exist('traj_type','var') == 0
    traj_type = 'infty';
end

traj_type = 'infty';

failure.type = 'all';
failure.amplitude = 0.1;
failure.amplitude_2 = 0.1;

bridge_type = 'cubic';

time_infor.val_length=250000;
val_length_rm = time_infor.val_length;

rmse_start_time = round(val_length_rm * 3/5);
rmse_end_time = time_infor.val_length - 100;

if strcmp(traj_type, 'circle') == 1
    frequency_set = round(linspace(10, 500, 15));
elseif strcmp(traj_type, 'infty') == 1
    frequency_set = round(linspace(10, 500, 15) / 2);
end

iteration = 50;

failure.type = 'none';

rmse_set = zeros(length(frequency_set), iteration);

idx=1;

for f_idx = 1:length(frequency_set)
    rmse_parfor_set = zeros(1, iteration);
    for repeat_i = 1:iteration
        save_rend=0;
        load('./save_file/all_traj_06282022.mat')
        time_infor.val_length = val_length_rm;
        traj_frequency = frequency_set(f_idx);

        val_and_update;

        rmse_parfor_set(repeat_i) = func_rmse(data_pred, data_control, rmse_start_time, rmse_end_time);
        aaa = 1;
    end
    rmse_parfor_set = sort(rmse_parfor_set);
    rmse_set(f_idx, :) = rmse_parfor_set;
end


save_speed_iter.traj_type = traj_type;

save_speed_iter.(['frequency_set', num2str(idx)]) = frequency_set;
save_speed_iter.(['rmse_set', num2str(idx)]) = rmse_set;

save(['save_data/save_speed_iter_' traj_type, '_' time_today, '_' num2str(randi(999)) '.mat'], "save_speed_iter")
































