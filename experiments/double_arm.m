%% clear all

disp('Preparing for clearing all and restart...')
m=input('Do you want to continue, Y/N [Y]:', 's');
if m ~= 'Y'
    return
end

clear all
close all
clc

addpath('./tools/')

%% pre training
dt=0.01;
% input_infor={'q'};
input_infor={'xy', 'qdt'};
% input_infor={'q', 'qdt'};
% input_infor={'xy', 'qdt', 'q2dt'};

% in-out dimension
dim_in=length(input_infor) * 4;
dim_out=2;

% double robot arm properties
m1=1;m2=1;
l1=0.5;l2=0.5;
lc1=0.25;lc2=0.25;
I1=0.03;I2=0.03;
properties=[m1, m2, l1, l2, lc1, lc2, I1, I2];

reset_t=80; 
train_t = 200000;
val_t=500;
noise_level=2.0*10^(-2);
% add pertubations and measurement noise to the model(and uncertainty)
% proportional to the real value.
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

%% training
tic;
[Wout, r_end] = func_reservoir_train(data_reservoir, time_infor, input_infor, res_infor, dim_in, dim_out);
toc;

clearvars data_reservoir
%% load data
load_data = 1;
if load_data==1
    load('./save_file/all_traj_06282022.mat')
end

%% validating

rng('shuffle')
% write the code that can update and pause
% traj_type = 'infty';
% traj_type = 'lorenz';
% traj_type = 'mg17';
% traj_type = 'mg30';
% traj_type = 'circle'

failure.type = 'none';
blur.blur = 0;

plot_val_and_update=1;

disturbance = 0.05;
measurement_noise = 0.05;
plot_movie = 0;

traj_type = 'lorenz';
bridge_type = 'cubic';

time_infor.val_length=250000;

save_rend=0;
idx=1;

val_and_update;

% idx = 2
traj_type = 'circle';
bridge_type = 'cubic';
traj_frequency = 150;

time_infor.val_length=250000;
save_rend=1;
idx=2;

val_and_update;

% idx = 3
traj_type = 'mg17';
bridge_type = 'cubic';

time_infor.val_length=250000;
save_rend=1;
idx=3;

val_and_update;

% idx = 4
traj_type = 'infty';
bridge_type = 'cubic';
traj_frequency = 100;

time_infor.val_length=250000;
save_rend=1;
idx=4;

val_and_update;

% normal method to calculate rmse
rmse_length = 1000/dt;
error = abs(data_control(1:rmse_length, :) - data_pred(1:rmse_length, :));
rmse = sqrt(mean(mean(error.^2, 2)));


%% plot figure

% start_time = 1;
% end_time = 290000;
% 
% % plot trajectory
% figure();
% hold on
% plot(data_control(start_time:end_time, 1), data_control(start_time:end_time, 2),'r');
% plot(data_pred(start_time:end_time, 1), data_pred(start_time:end_time, 2),'b--');
% xlabel('x')
% ylabel('y')
% line([0, 0], [-1, 1], 'Color', 'black', 'LineStyle', '--')
% line([-1, 1], [0, 0], 'Color', 'black', 'LineStyle', '--')
% xlim([-1, 1])
% ylim([-1, 1])
% legend('desired trajectory', 'pred trajectory')
% % 
% % figure();
% % hold on
% % plot(data_control(100000:120000, 1), data_control(100000:120000, 2),'r');
% % plot(data_pred(100000:120000, 1), data_pred(100000:120000, 2),'b--');
% % xlabel('x')
% % ylabel('y')
% % line([0, 0], [-1, 1], 'Color', 'black', 'LineStyle', '--')
% % line([-1, 1], [0, 0], 'Color', 'black', 'LineStyle', '--')
% % xlim([-1, 1])
% % ylim([-1, 1])
% % legend('desired trajectory', 'pred trajectory')
% 
% % plot q
% q_control_plot=mod(q_control, pi);
% q_pred_plot=mod(q_pred, pi);
% 
% figure();
% hold on
% plot(q_control_plot(start_time:end_time,1), 'r')
% plot(q_pred_plot(start_time:end_time,1), 'b')
% xlabel('time step')
% ylabel('q(1)')
% legend('desired', 'pred')
% 
% figure()
% hold on
% plot(q_control_plot(start_time:end_time,2), 'r')
% plot(q_pred_plot(start_time:end_time,2), 'b')
% xlabel('time step')
% ylabel('q(2)')
% legend('desired', 'pred')
% 
% % plot dq/dt
% figure()
% hold on
% plot(qdt_control(start_time:end_time,1), 'r')
% plot(qdt_pred(start_time:end_time,1), 'b')
% xlabel('time step')
% ylabel('dq/dt(1)')
% legend('desired', 'pred')
% 
% figure()
% hold on
% plot(qdt_control(start_time:end_time,2), 'r')
% plot(qdt_pred(start_time:end_time,2), 'b')
% xlabel('time step')
% ylabel('dq/dt(2)')
% legend('desired', 'pred')
% 
% % plot d2q/dt2
% remove_transient = 1000;
% 
% figure()
% hold on
% plot(q2dt_control(start_time+remove_transient:end_time,1), 'r')
% plot(q2dt_pred(start_time+remove_transient:end_time,1), 'b')
% xlabel('time step')
% ylabel('d2q/dt2(1)')
% legend('desired', 'pred')
% 
% figure()
% hold on
% plot(q2dt_control(start_time+remove_transient:end_time,2), 'r')
% plot(q2dt_pred(start_time+remove_transient:end_time,2), 'b')
% xlabel('time step')
% ylabel('d2q/dt2(2)')
% legend('desired', 'pred')
% 
% figure()
% hold on
% plot(tau_control(start_time+remove_transient:end_time,1), 'r')
% plot(tau_pred(start_time+remove_transient:end_time,1), 'b')
% xlabel('time step')
% ylabel('tau(1)')
% legend('desired', 'pred')
% 
% figure()
% hold on
% plot(tau_control(start_time+remove_transient:end_time,2), 'r')
% plot(tau_pred(start_time+remove_transient:end_time,2), 'b')
% xlabel('time step')
% ylabel('tau(2)')
% legend('desired', 'pred')
% 
% % figure()
% % hold on
% % plot(tau_control(start_time:end_time,1), 'r')
% % plot(tau_control(start_time:end_time,2), 'b')
% % xlabel('time step')
% % ylabel('tau')
% % legend('tau1', 'tau2')
% 
% % plot movie for validation
% plot_movie_val = 0;
% start_step=start_time;
% movie_step=500;
% time_all=end_time;
% line_property='dotted';
% q1=q_pred(:, 1);
% q2=q_pred(:, 2);
% if plot_movie_val == 1
%     func_plot_movie(start_step, movie_step, time_all, q1, q2, properties, line_property)
% end


%% save data

% save('./save_file/8d_input_04222022_lorenz_recover_from_failure.mat')
% save('./save_file/8d_input_04222022_mg17.mat')
% save('./save_file/reset_to_same_point_05182022_infty.mat')
% save('./save_file/reset_to_random_point_05242022_infty.mat')

% save('./save_file/reto_random_work_for_all_traj_05242022.mat')
% 
% save('./save_data/4_trajectory_08182022.mat', "save_all_traj")

% save('./save_file/all_traj_06282022.mat', 'time_infor', 'input_infor', 'res_infor', 'properties', 'dim_in', 'dim_out', 'Wout', 'r_end', 'dt', 'reset_t', 'noise_level')
% save('./save_file/all_traj_06282022_SI.mat', 'data_reservoir')


%% success test

% aa = 1;
% 
% plot_movie = 0;
% traj_type = 'infty';
% bridge_type = 'cubic';
% plot_val_and_update=0;
% 
% time_infor.val_length=150000;
% rmse_start_time = round(time_infor.val_length * 2/5);
% rmse_end_time = time_infor.val_length - 100;
% 
% iteration = 50;
% rmse_parfor_set = zeros(1, iteration);
% for repeat_i = 1:iteration
%     load('./save_file/all_traj_06282022.mat')
%     val_and_update;
%     data_pred=output_infor.data_pred;
%     data_control=control_infor.data_control;
%     rmse_parfor_set(repeat_i) = func_rmse(data_pred, data_control, rmse_start_time, rmse_end_time);
% end
% 
% rmse_parfor_set = sort(rmse_parfor_set);




















