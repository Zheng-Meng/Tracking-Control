% validate
if save_rend == 0
    start_info.q=0;
    start_info.qdt=0;
    start_info.q2dt=0;
    start_info.tau=0;

    r_end = zeros(size(res_infor.res_net, 1), 1);
end

if exist('traj_frequency','var') == 0
    if strcmp(traj_type, 'lorenz') == 1
        traj_frequency = 100;
    elseif strcmp(traj_type, 'cirlce') == 1
        traj_frequency = 150;
    else
        traj_frequency = 75;
    end
end

rng('shuffle')

[control_infor, output_infor, time_infor, r_end] = func_reservoir_validate(traj_type,...
bridge_type, time_infor, input_infor, res_infor, start_info, properties, dim_in, dim_out, ...
    Wout, r_end, dt, plot_movie, save_rend,failure,blur,traj_frequency);

% update
data_pred=output_infor.data_pred;
q_pred=output_infor.q_pred;
qdt_pred=output_infor.qdt_pred;
q2dt_pred=output_infor.q2dt_pred;
tau_pred=output_infor.tau_pred;

q_control=control_infor.q_control;
qdt_control=control_infor.qdt_control;
q2dt_control=control_infor.q2dt_control;
tau_control=control_infor.tau_control;
data_control=control_infor.data_control;

val_length=time_infor.val_length;

start_info.q=q_pred(val_length-3,:);
start_info.qdt=qdt_pred(val_length-3, :);
start_info.q2dt=q2dt_pred(val_length-3,:);
start_info.tau=tau_pred(val_length-3,:);

save_all_traj.(['control_', num2str(idx)]) = control_infor;
save_all_traj.(['output_', num2str(idx)]) = output_infor;

% plot trajectory

start_time=1;
end_time=val_length-100;

if exist('plot_val_and_update','var') == 0
    plot_val_and_update = 0;
end

if plot_val_and_update==1
    figure();
    hold on
    plot(data_control(start_time:end_time, 1), data_control(start_time:end_time, 2),'r');
    plot(data_pred(start_time:end_time, 1), data_pred(start_time:end_time, 2),'b--');
    xlabel('x')
    ylabel('y')
    line([0, 0], [-1, 1], 'Color', 'black', 'LineStyle', '--')
    line([-1, 1], [0, 0], 'Color', 'black', 'LineStyle', '--')
    xlim([-1, 1])
    ylim([-1, 1])
    legend('desired trajectory', 'pred trajectory')
end



















