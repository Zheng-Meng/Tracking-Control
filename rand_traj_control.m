clear all
close all
clc

addpath('./tools/')

% load training data
load('./save_file/all_traj_06282022.mat')

%% 

traj_set = ["infty", "circle", "astroid", "fermat", ...
    "lissajous", "talbot", "heart", "lorenz", "chua", "rossler", ...
    "sprott_1", "sprott_4", "mg17", "mg30", "epitrochoid"];

order = randperm(length(traj_set));
traj_set= traj_set(order);

plot_val_and_update = 1;
disturbance = 0.1;
measurement_noise = 0.1;
plot_movie = 0;
bridge_type = 'cubic';
failure.type = 'none';
blur.blur = 0;
idx=1;

val_length_all=150000;

for ii = 1:length(traj_set)
    rng('shuffle')
    traj_type = traj_set(ii);
    time_infor.val_length=val_length_all;

    if ii == 1
        save_rend=0;
    else
        save_rend=1;
    end

    if exist('traj_frequency','var') == 1
        if strcmp(traj_type, 'lorenz') == 1
            traj_frequency = 100;
        elseif strcmp(traj_type, 'cirlce') == 1
            traj_frequency = 150;
        else
            traj_frequency = 75;
        end
    end

    val_and_update;

    idx = idx+1;
    if ii == 1
        aaa = 1;
    end
end

time_infor.val_length=val_length_all;

time_today = datestr(now, 'mmddyyyy');
% save(['./save_data/15traj_', time_today, '_', num2str(randi(999)), '.mat'], "val_length_all", "save_all_traj", "traj_set")









































