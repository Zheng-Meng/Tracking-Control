function [] = func_double_arm()

rng('shuffle')

time_today = datestr(now, 'mmddyyyy');

dt=0.01;
input_infor={'xy', 'qdt'};

% in-out dimension
dim_in=length(input_infor) * 4;
dim_out=2;

% double robot arm properties
m1=1;m2=1;
l1=0.5;l2=0.5;
lc1=0.25;lc2=0.25;
I1=0.03;I2=0.03;
properties=[m1, m2, l1, l2, lc1, lc2, I1, I2];

reset_t=80; % 70 to get good results for infty symbol, and noise level:0.8*0.5*10^(-2)
train_t = 200000;
val_t=500;

noise_level=2.0*10^(-2);
disturbance = 0.00;
measurement_noise = 0.00;

n=200;
hyperpara_set = [0.756250, 0.756250, 0.843750, -3.125, 106.71875, 2.0];
eig_rho = hyperpara_set(1);
W_in_a = hyperpara_set(2);
alpha = hyperpara_set(3);
beta = 10^hyperpara_set(4);
k = round( hyperpara_set(5)/200*n);
kb = hyperpara_set(6);

W_in = W_in_a*(2*rand(n,dim_in)-1);
res_net=sprandsym(n,k/n);
eig_D=eigs(res_net,1);
res_net=(eig_rho/(abs(eig_D))).*res_net;
res_net=full(res_net);

section_len=round(reset_t/dt); % 30
washup_length=round(1002/dt);
train_length=round(train_t/dt)+round(5/dt); % 50000
val_length=round(val_t/dt); % 300
time_length = train_length + 2 * val_length + 3 * washup_length + 100;

res_infor=struct('W_in', W_in, 'res_net', res_net, 'alpha', alpha, 'kb', kb, 'beta', beta, 'n', n);
time_infor=struct('section_len', section_len, 'washup_length', washup_length, ...
    'train_length', train_length, 'val_length', val_length, 'time_length', time_length);

% generate training and validation data
[xy, q, qdt, q2dt, tau] = robot_data_generator(time_infor, noise_level, dt, properties);
xy=xy(washup_length:end, :);
q=q(washup_length:end, :);
qdt=qdt(washup_length:end, :);
q2dt=q2dt(washup_length:end, :);
tau=tau(washup_length:end, :);
data_reservoir = struct('xy', xy, 'q', q, 'qdt', qdt, 'q2dt', q2dt, 'tau', tau);

clearvars xy q qdt q2dt tau

tic;
[Wout, r_end] = func_reservoir_train(data_reservoir, time_infor, input_infor, res_infor, dim_in, dim_out);
toc;

clearvars data_reservoir

rng('shuffle')

failure.type = 'none';
blur.blur = 0;

plot_val_and_update=0;

disturbance = 0.00;
measurement_noise = 0.00;

rmse_start_time = 200000;
rmse_end_time = 300000-100;

plot_movie = 0;
traj_type = 'lorenz';
bridge_type = 'cubic';

time_infor.val_length=300000;

save_rend=0;
idx=1;

val_and_update;

rmse_1 = func_rmse(data_pred, data_control, rmse_start_time, rmse_end_time);

% idx = 2
traj_type = 'circle';
bridge_type = 'cubic';

time_infor.val_length=300000;
save_rend=1;
idx=2;

val_and_update;

rmse_2 = func_rmse(data_pred, data_control, rmse_start_time, rmse_end_time);

% idx = 3
traj_type = 'mg17';
bridge_type = 'cubic';

time_infor.val_length=300000;
save_rend=1;
idx=3;

val_and_update;

rmse_3 = func_rmse(data_pred, data_control, rmse_start_time, rmse_end_time);

% idx = 41
traj_type = 'infty';
bridge_type = 'cubic';

time_infor.val_length=300000;
save_rend=1;
idx=4;

val_and_update;

rmse_4 = func_rmse(data_pred, data_control, rmse_start_time, rmse_end_time);

a_rmse = (rmse_1 + rmse_2 + rmse_3 + rmse_4) / 4;

save(['./choose_file/all_traj_', time_today, '_' num2str(randi(9999)) '_' num2str(randi(9999)) '.mat'], 'time_infor', 'input_infor', 'res_infor', 'properties', 'dim_in', 'dim_out', 'Wout', 'r_end', 'dt', 'reset_t', 'noise_level', 'a_rmse')

end








