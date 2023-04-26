clear all
close all
clc

addpath('./tools/')

load('./save_file/all_traj_06282022.mat')

% choose the reference trajectory
traj_type = 'lorenz';
% traj_type = 'circle';
% traj_type = 'mg17';
% traj_type = 'infty';

% traj_frequency = 75;

time_infor.val_length=150000;
bridge_type = 'cubic';
failure.type = 'none';
disturbance = 0.00;
measurement_noise = 0.00;
plot_movie = 0;
blur.blur = 0;

save_rend = 0;
idx = 0;

plot_val_and_update = 1;

% let the well-trained machine to follow the given reference.
val_and_update



















