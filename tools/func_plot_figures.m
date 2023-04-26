function [] = func_plot_figures(control_infor, output_infor, val_length, plot_movie)

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

start_time = 1;
end_time = val_length - 10000;

% plot trajectory
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

% plot q
q_control_plot=mod(q_control, pi);
q_pred_plot=mod(q_pred, pi);

figure();
hold on
plot(q_control_plot(start_time:end_time,1), 'r')
plot(q_pred_plot(start_time:end_time,1), 'b')
xlabel('time step')
ylabel('q(1)')
legend('desired', 'pred')

figure()
hold on
plot(q_control_plot(start_time:end_time,2), 'r')
plot(q_pred_plot(start_time:end_time,2), 'b')
xlabel('time step')
ylabel('q(2)')
legend('desired', 'pred')

% plot dq/dt
figure()
hold on
plot(qdt_control(start_time:end_time,1), 'r')
plot(qdt_pred(start_time:end_time,1), 'b')
xlabel('time step')
ylabel('dq/dt(1)')
legend('desired', 'pred')

figure()
hold on
plot(qdt_control(start_time:end_time,2), 'r')
plot(qdt_pred(start_time:end_time,2), 'b')
xlabel('time step')
ylabel('dq/dt(2)')
legend('desired', 'pred')

% plot d2q/dt2
remove_transient = 1000;

figure()
hold on
% plot(q2dt_control(start_time+remove_transient:end_time,1), 'r')
plot(q2dt_pred(start_time+remove_transient:end_time,1), 'b')
xlabel('time step')
ylabel('d2q/dt2(1)')
legend('pred')

figure()
hold on
% plot(q2dt_control(start_time+remove_transient:end_time,2), 'r')
plot(q2dt_pred(start_time+remove_transient:end_time,2), 'b')
xlabel('time step')
ylabel('d2q/dt2(2)')
legend('pred')

figure()
hold on
% plot(tau_control(start_time+remove_transient:end_time,1), 'r')
plot(tau_pred(start_time+remove_transient:end_time,1), 'b')
xlabel('time step')
ylabel('tau(1)')
legend('pred')

figure()
hold on
% plot(tau_control(start_time+remove_transient:end_time,2), 'r')
plot(tau_pred(start_time+remove_transient:end_time,2), 'b')
xlabel('time step')
ylabel('tau(2)')
legend('pred')

% figure()
% hold on
% plot(tau_control(start_time:end_time,1), 'r')
% plot(tau_control(start_time:end_time,2), 'b')
% xlabel('time step')
% ylabel('tau')
% legend('tau1', 'tau2')

% plot movie for validation
plot_movie_val = plot_movie;
start_step=start_time;
movie_step=500;
time_all=end_time;
line_property='dotted';
q1=q_pred(:, 1);
q2=q_pred(:, 2);
if plot_movie_val == 1
    func_plot_movie(start_step, movie_step, time_all, q1, q2, properties, line_property)
end




end