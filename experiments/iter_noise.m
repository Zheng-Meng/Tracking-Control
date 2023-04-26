%
% disturbance = 1.00;
% measurement_noise = 0.00;

% try for only the first one and then only the second one
% then try them together
% then try uncertainty

time_today = datestr(now, 'mmddyyyy');

%% only disturbance
measurement_noise = 0.00;

% disturbance_set = exp(linspace(log(0.01), log(10), 20));
disturbance_set = linspace(1, 5, 20);
error_set = zeros(1, length(disturbance_set));

rmse_length = 1000/dt;
weight=zeros(rmse_length,1);
% alpha_w=0.2;
for ii=1:rmse_length
    weight(ii)=ii;
end
weight = weight/norm(weight, 1);

for di_id = 1:length(disturbance_set)
    disturbance = disturbance_set(di_id);
    plot_movie = 0;
    traj_type = 'infty';
    % traj_type = 'lorenz';
    % traj_type = 'mg17';
    % traj_type = 'mg30';
    bridge_type = 'cubic';

    time_infor.val_length=120000;

    [control_infor, output_infor, time_infor] = func_reservoir_validate(traj_type,...
        bridge_type, time_infor, input_infor, res_infor, properties, dim_in, dim_out, ...
        Wout, dt, plot_movie, disturbance, measurement_noise);

    data_pred=output_infor.data_pred;
    q_pred=output_infor.q_pred;
    qdt_pred=output_infor.qdt_pred;
    q2dt_pred=output_infor.q2dt_pred;
    tau_pred=output_infor.tau_pred;

    q_control=control_infor.q_control;
    qdt_control=control_infor.qdt_control;
    q2dt_control=control_infor.q2dt_control;
    tau_control=control_infor.tau_control;
    % q_control_all=control_infor.q_control_all;
    % qdt_control_all=control_infor.qdt_control_all;
    data_control=control_infor.data_control;

    val_length=time_infor.val_length;

    % calculate rmse
    error = abs(data_control(1:rmse_length, :) - data_pred(1:rmse_length, :));
    error = sum(weight.*mean(error.^2, 2));

    error_set(di_id) = error;
end

nan_id = isnan(error_set);
error_set_plot = error_set;
error_set_plot(nan_id) = 5*max(error_set, [], 'omitnan');

save_data.disturbance_set = disturbance_set;
save_data.disturbance_error_set = error_set;
save_data.disturbance_error_set_plot = error_set_plot;

save(['./save_data/disturbance_result_normscale_', time_today, '.mat'], "save_data")


%% for plot

figure();
% hold on
% semilogx(disturbance_set, error_set_plot, 'o-')
plot(disturbance_set, error_set_plot, 'o-')
xlabel('gaussian disturbance, \sigma')
ylabel('error')


%% only measurement
disturbance = 0.00;
measurement_set = exp(linspace(log(0.01), log(10), 20));
% measurement_set = linspace(1, 5, 20);
error_set = zeros(1, length(measurement_set));

rmse_length = 1000/dt;
weight=zeros(rmse_length,1);
for ii=1:rmse_length
    weight(ii)=ii;
end
weight = weight/norm(weight, 1);

for di_id = 1:length(measurement_set)
    measurement_noise = measurement_set(di_id);
    plot_movie = 0;
    traj_type = 'infty';
    % traj_type = 'lorenz';
    % traj_type = 'mg17';
    % traj_type = 'mg30';
    bridge_type = 'cubic';

    time_infor.val_length=120000;

    [control_infor, output_infor, time_infor] = func_reservoir_validate(traj_type,...
        bridge_type, time_infor, input_infor, res_infor, properties, dim_in, dim_out, ...
        Wout, dt, plot_movie, disturbance, measurement_noise);

    data_pred=output_infor.data_pred;
    q_pred=output_infor.q_pred;
    qdt_pred=output_infor.qdt_pred;
    q2dt_pred=output_infor.q2dt_pred;
    tau_pred=output_infor.tau_pred;

    q_control=control_infor.q_control;
    qdt_control=control_infor.qdt_control;
    q2dt_control=control_infor.q2dt_control;
    tau_control=control_infor.tau_control;
    % q_control_all=control_infor.q_control_all;
    % qdt_control_all=control_infor.qdt_control_all;
    data_control=control_infor.data_control;

    val_length=time_infor.val_length;

    % calculate rmse
    error = abs(data_control(1:rmse_length, :) - data_pred(1:rmse_length, :));
    error = sum(weight.*mean(error.^2, 2));

    error_set(di_id) = error;
end

nan_id = isnan(error_set);
error_set_plot = error_set;
error_set_plot(nan_id) = 5*max(error_set, [], 'omitnan');

save_data.measurement_set = measurement_set;
save_data.measurement_error_set = error_set;
save_data.measurement_error_set_plot = error_set_plot;

save(['./save_data/measurement_result_logscale_', time_today, '.mat'], "save_data")

%% for plot

load('./save_data/measurement_result_logscale_05202022.mat')

error_set_plot = save_data.measurement_error_set;

figure();
% hold on
semilogx(measurement_set, error_set_plot, 'o-')
% plot(measurement_set, error_set_plot, 'o-')
xlabel('gaussian measurement noise, \sigma')
ylabel('error')



%% disturbance and measurement noise

disturbance_set = exp(linspace(log(0.01), log(10), 20));
measurement_set = exp(linspace(log(0.01), log(10), 20));

error_set = zeros(length(disturbance_set), length(measurement_set));

rmse_length = 1000/dt;
weight=zeros(rmse_length,1);
for ii=1:rmse_length
    weight(ii)=ii;
end
weight = weight/norm(weight, 1);

for di_id = 1:length(disturbance_set)
    for mn_id = 1:length(measurement_set)
        disturbance = disturbance_set(di_id);
        measurement_noise = measurement_set(mn_id);
        plot_movie = 0;
        traj_type = 'infty';
        % traj_type = 'lorenz';
        % traj_type = 'mg17';
        % traj_type = 'mg30';
        bridge_type = 'cubic';

        time_infor.val_length=120000;

        [control_infor, output_infor, time_infor] = func_reservoir_validate(traj_type,...
            bridge_type, time_infor, input_infor, res_infor, properties, dim_in, dim_out, ...
            Wout, dt, plot_movie, disturbance, measurement_noise);

        data_pred=output_infor.data_pred;
        q_pred=output_infor.q_pred;
        qdt_pred=output_infor.qdt_pred;
        q2dt_pred=output_infor.q2dt_pred;
        tau_pred=output_infor.tau_pred;

        q_control=control_infor.q_control;
        qdt_control=control_infor.qdt_control;
        q2dt_control=control_infor.q2dt_control;
        tau_control=control_infor.tau_control;
        % q_control_all=control_infor.q_control_all;
        % qdt_control_all=control_infor.qdt_control_all;
        data_control=control_infor.data_control;

        val_length=time_infor.val_length;

        % calculate rmse
        error = abs(data_control(1:rmse_length, :) - data_pred(1:rmse_length, :));
        error = sum(weight.*mean(error.^2, 2));

        error_set(di_id, mn_id) = error;
    end
end

nan_id = isnan(error_set);
error_set_plot = error_set;
error_set_plot(nan_id) = 3*max(max(error_set));

save_data.disturbance_set_heat = disturbance_set;
save_data.measurement_set_heat = measurement_set;
save_data.heat_error_set = error_set;
save_data.heat_error_set_plot = error_set_plot;

save(['./save_data/heat_result_logscale_', time_today, '.mat'], "save_data")

%% for plot

figure();
surf(disturbance_set, measurement_set, error_set_plot);
xlabel('measurement noise, \sigma')
ylabel('disturbance noise, \sigma')
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
colorbar
view(0, 90)
title('error')



















