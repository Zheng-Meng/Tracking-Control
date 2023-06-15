function [control_infor, time_infor] = func_desired_traj(traj_type, bridge_type, time_infor, control_infor, properties, dt, plot_movie, traj_frequency)

plot_figures_inside = 0; 

val_length=time_infor.val_length;
train_length=time_infor.train_length;
% Get the two-arms property
[m1, m2, l1, l2, lc1, lc2, I1, I2] = matsplit(properties);
% Get the initial information
q_control=control_infor.q_control;
qdt_control=control_infor.qdt_control;
q2dt_control=control_infor.q2dt_control;
tau_control=control_infor.tau_control;

value_q=q_control(1,:);
x_start=l1*cos(value_q(1))+l2*cos(value_q(1)+value_q(2));
y_start=l1*sin(value_q(1))+l2*sin(value_q(1)+value_q(2));

t=0:dt:val_length;

%% Different traj_types, the readers can add their own trajectories
if strcmp(traj_type, 'infty')==1
    x = 0.25*sin(2*pi*t*1/(2*traj_frequency));
    y = 0.15*sin(2*pi*t*1/traj_frequency);  

    x = x(1:2*val_length+1);
    y = y(1:2*val_length+1);
elseif strcmp(traj_type, 'circle')
    x = 0.5 * cos(2 * pi * t * 1 / traj_frequency);
    y = 0.5 * sin(2 * pi * t * 1 / traj_frequency);

    x = x(1:2*val_length+1);
    y = y(1:2*val_length+1);
elseif strcmp(traj_type, 'astroid')
    x = 0.4 * (cos(2 * pi * t * 1 / 250).^3);
    y = 0.4 * (sin(2 * pi * t * 1 / 250).^3);
    
    x = x(1:2*val_length+1);
    y = y(1:2*val_length+1);

elseif strcmp(traj_type, 'heart')
    f = 250;
    a = 0.4;
    theta = 2 * pi * t * 1 / f;

    r = 1 - sin(theta);

    x = r.*cos(theta);
    y = r.*sin(theta);

    x = normalize(x, 'range', [-a, a]);
    y = normalize(y, 'range', [-a, a]);

    x = x(1:2*val_length+1);
    y = y(1:2*val_length+1);

elseif strcmp(traj_type, 'epitrochoid')
    a = 5;
    b = 3;
    c = 5;

    f = 200;
    
    x = (a + b) * cos(2 * pi * t * 1 / f) - c * cos((a/b + 1) * 2 * pi * t * 1 / f);
    y = (a + b) * sin(2 * pi * t * 1 / f) - c * sin((a/b + 1) * 2 * pi * t * 1 / f);

    x = normalize(x, 'range', [-0.4, 0.4]);
    y = normalize(y, 'range', [-0.4, 0.4]);

    x = x(1:2*val_length+1);
    y = y(1:2*val_length+1);

elseif strcmp(traj_type, 'fermat')
    f = 100;
    a = 0.5;
    theta = 2 * pi * t * 1 / f;

    r = sqrt(a^2 * theta);

    x = r.*cos(theta);
    y = r.*sin(theta);

    x = normalize(x, 'range', [-4, 4]);
    y = normalize(y, 'range', [-4, 4]);

    x = x(1:2*val_length+1);
    y = y(1:2*val_length+1);

elseif strcmp(traj_type, 'lissajous')
    f = 300;
    a = 1;
    b = 3;
    A = 1;
    B = 1;

    x = A * sin(a * t * 2 * pi / f  + pi / 4);
    y = B * sin(b * t * 2 * pi / f) ;

    xx = normalize(x, 'range', [-0.3, 0.3]);
    yy = normalize(y, 'range', [-0.3, 0.3]);
    
    x = yy(1:2*val_length+1);
    y = xx(1:2*val_length+1);

elseif strcmp(traj_type, 'talbot')
    f = 400;
    a = 1.1;
    b = 0.666;
    ff = 1;

    x = (a^2 + ff^2 .* sin(t * 2 * pi / f).^2) .* cos(t * 2 * pi / f) / a;
    y = (a^2 - 2*ff^2 + ff^2 .* sin(t * 2 * pi / f).^2) .* sin(t * 2 * pi / f) / b;

    x = normalize(x, 'range', [-0.4, 0.4]);
    y = normalize(y, 'range', [-0.3, 0.3]);

    x = x(1:2*val_length+1);
    y = y(1:2*val_length+1);

elseif strcmp(traj_type, 'lorenz')
    load('./read_data/lorenz.mat')
    lorenz_xy=ts_train(1000:1000+round(val_length*2.1), 1:2);

    lorenz_x = normalize(lorenz_xy(:, 1), 'range', [-0.5, 0.5]);
    lorenz_y = normalize(lorenz_xy(:, 2), 'range', [-0.5, 0.5]);

    x = lorenz_x(1:val_length*2+1, :)';
    y = lorenz_y(1:val_length*2+1, :)';

    y = y - 0.3;

elseif strcmp(traj_type, 'chua')
    load('./read_data/chua.mat')
    chua_xy = ts_train(10000:10000+round(val_length*2.1), 1:2);

    chua_x = normalize(chua_xy(:, 1), 'range', [-0.5, 0.5]);
    chua_y = normalize(chua_xy(:, 2), 'range', [-0.3, 0.3]);

    chua_x = interp1(1:length(chua_x), chua_x, 0:0.3:length(chua_x), 'pchip');
    chua_y = interp1(1:length(chua_y), chua_y, 0:0.3:length(chua_y), 'pchip');

    x = chua_x(1:val_length*2+1)';
    y = chua_y(1:val_length*2+1)';

    x = reshape(x, 1, []);
    y = reshape(y, 1, []);

elseif strcmp(traj_type, 'rossler')
    load('./read_data/rossler.mat')
    rossler_xy = ts_train(10000:10000+round(val_length*2.1), 1:2);

    rossler_x = normalize(rossler_xy(:, 1), 'range', [-0.35, 0.35]);
    rossler_y = normalize(rossler_xy(:, 2), 'range', [-0.35, 0.35]);

    rossler_x = interp1(1:length(rossler_x), rossler_x, 0:0.55:length(rossler_x), 'pchip');
    rossler_y = interp1(1:length(rossler_y), rossler_y, 0:0.55:length(rossler_y), 'pchip');

    x = rossler_x(1:val_length*2+1)';
    y = rossler_y(1:val_length*2+1)';

    x = reshape(x, 1, []);
    y = reshape(y, 1, []);

elseif strcmp(traj_type, 'sprott_1')
    load('./read_data/sprott_1.mat')
    sprott_xy = ts_train(10000:10000+round(val_length*2.1), 1:2);

    sprott_x = normalize(sprott_xy(:, 1), 'range', [-0.4, 0.4]);
    sprott_y = normalize(sprott_xy(:, 2), 'range', [-0.4, 0.4]);

    sprott_x = interp1(1:length(sprott_x), sprott_x, 0:0.6:length(sprott_x), 'pchip');
    sprott_y = interp1(1:length(sprott_y), sprott_y, 0:0.6:length(sprott_y), 'pchip');

    x = sprott_x(1:val_length*2+1)';
    y = sprott_y(1:val_length*2+1)';

    x = reshape(x, 1, []);
    y = reshape(y, 1, []);

elseif strcmp(traj_type, 'sprott_4')
    load('./read_data/sprott_4.mat')
    sprott_xy = ts_train(10000:10000+round(val_length*2.1), 1:2);

    sprott_x = normalize(sprott_xy(:, 1), 'range', [-0.4, 0.4]);
    sprott_y = normalize(sprott_xy(:, 2), 'range', [-0.4, 0.4]);

    sprott_x = interp1(1:length(sprott_x), sprott_x, 0:0.6:length(sprott_x), 'pchip');
    sprott_y = interp1(1:length(sprott_y), sprott_y, 0:0.6:length(sprott_y), 'pchip');

    x = sprott_x(1:val_length*2+1)';
    y = sprott_y(1:val_length*2+1)';

    x = reshape(x, 1, []);
    y = reshape(y, 1, []);


elseif strcmp(traj_type, 'mg17')
    load('./read_data/MG17.mat')
    mg=ts_train(10000:10000+round(val_length*2.1));

    mg_x = normalize(mg(10000-1700:end-1700), 'range', [-0.5, 0.5]);
    mg_y = normalize(mg(10000:end), 'range', [-0.5, 0.5]);

    mg_x = interp1(1:length(mg_x), mg_x, 0:0.25:length(mg_x), 'pchip');
    mg_y = interp1(1:length(mg_y), mg_y, 0:0.25:length(mg_y), 'pchip');

    x = mg_x(1:val_length*2+1)';
    y = mg_y(1:val_length*2+1)';

    x = reshape(x, 1, []);
    y = reshape(y, 1, []);

    x = x - 0.1;
    y = y - 0.1;
elseif strcmp(traj_type, 'mg30')
    load('./read_data/MG30.mat')
    mg=ts_train(10000:10000+round(val_length*2.1));
    
    mg_x = normalize(mg(10000:end), 'range', [-0.4, 0.4]);
    mg_y = normalize(mg(10000-3000:end-3000), 'range', [-0.4, 0.4]);
   
    mg_x = interp1(1:length(mg_x), mg_x, 0:0.25:length(mg_x), 'pchip');
    mg_y = interp1(1:length(mg_y), mg_y, 0:0.25:length(mg_y), 'pchip');
    
    x = mg_x(1:val_length*2+1)';
    y = mg_y(1:val_length*2+1)';
    
    x = reshape(x, 1, []);
    y = reshape(y, 1, []);

    x = x + 0.25;
    y = y;

elseif strcmp(traj_type, 'lorenz96')
    load('./read_data/lorenz96.mat')

    lorenz_highD = ts_train(10000:10000+round(val_length*2.1), :);

    lorenz_x = normalize(lorenz_highD(:, 1), 'range', [-0.6, 0.6]);
    lorenz_y = normalize(lorenz_highD(:, 2), 'range', [-0.7, 0.7]);

    x = lorenz_x(1:val_length*2+1)';
    y = lorenz_y(1:val_length*2+1)';

    x = reshape(x, 1, []);
    y = reshape(y, 1, []);

    x = x + 0.1;
    y = y;
else
    disp('error: please input the traj type!')
end

%% Build bridge to the reference trajectory
add_id=1;
closet_value = Inf;
for i =1:min(100000, val_length)
    dis_value = sqrt((x_start-x(i)).^2+(y_start-y(i)).^2);
    if dis_value < closet_value
        if round(x(i), 6)==0 && round(y(i), 6)==0
            continue
        end
        closet_value = dis_value;
        add_id=i;
    end
end
bridge_point = [x(add_id), y(add_id)];
bridge_len = sqrt((bridge_point(1)-x_start).^2 + (bridge_point(2)-y_start).^2);
bridge_time = round(bridge_len * 1 / dt);
% We provide two methods to build the bridge: cubic and linear
if strcmp(bridge_type, 'cubic') == 1
    t_bg = 0:dt:bridge_time;
    theta0=q_control(1,:);
    theta0_dot=qdt_control(1,:);
    theta0=mod(theta0, 2*pi);
    if theta0(2) > pi
        theta0(2) = theta0(2)-2*pi;
    end
    q2_bg(1)=acos((x(add_id).^2+y(add_id).^2-l1^2-l2^2)/(2*l1*l2));
    q2_bg(2)=acos((x(add_id+1).^2+y(add_id+1).^2-l1^2-l2^2)/(2*l1*l2));
    if theta0(2)<0
        q2_bg=-q2_bg;
    end
    
    q1_bg(1)=atan(y(add_id)./x(add_id))-atan(l2*sin(q2_bg(1))./(l1+l2*cos(q2_bg(1))));
    q1_bg(2)=atan(y(add_id+1)./x(add_id+1))-atan(l2*sin(q2_bg(2))./(l1+l2*cos(q2_bg(2))));
    
    for ij = 1:length(q1_bg)
        if x(add_id + ij-1) < 0 && y(add_id + ij-1) > 0
            q1_bg(ij)=q1_bg(ij)+pi;
        elseif x(add_id + ij-1) < 0 && y(add_id + ij-1) < 0
            q1_bg(ij)=q1_bg(ij)+pi;
        elseif x(add_id + ij-1) > 0 && y(add_id + ij-1) < 0
            q1_bg(ij)=q1_bg(ij)+2*pi;
        end
    end
    
    theta1=[q1_bg(1), q2_bg(1)];
    
    x_back_test = l1*cos(theta1(1)) + l2*cos(theta1(1)+theta1(2));
    y_back_test = l1*sin(theta1(1)) + l2*sin(theta1(1)+theta1(2));
    
    theta1_dot=[(q1_bg(2)-q1_bg(1))/dt, (q2_bg(2)-q2_bg(1))/dt];
    
    a0=theta0;
    a1=theta0_dot;
    a2=3.*(theta1-theta0)./(bridge_time^2)-2.*theta0_dot./bridge_time-theta1_dot/bridge_time;
    a3=-2.*(theta1-theta0)./(bridge_time^3)+(theta1_dot+theta0_dot)/(bridge_time^2);
    
    q1_bridge=a0(1) + a1(1).*t_bg + a2(1).*(t_bg.^2) + a3(1).*(t_bg.^3);
    q2_bridge=a0(2) + a1(2).*t_bg + a2(2).*(t_bg.^2) + a3(2).*(t_bg.^3);
    
    truePosition=zeros(2, length(q1_bridge));
    truePosition(1,:)=l1*cos(q1_bridge)+l2*cos(q1_bridge+q2_bridge);
    truePosition(2,:)=l1*sin(q1_bridge)+l2*sin(q1_bridge+q2_bridge);
    
    x=[truePosition(1, 2:end-1), reshape(x(add_id:end), 1, [])];
    y=[truePosition(2, 2:end-1), reshape(y(add_id:end), 1, [])];
elseif strcmp(bridge_type, 'linear') == 1
    bridge = waypointTrajectory([x_start, y_start, 0; bridge_point(1), bridge_point(2), 0], 'TimeOfArrival', [0, bridge_time], 'SampleRate', round(bridge_time/dt));
    
    truePosition = zeros(bridge.SampleRate * bridge.TimeOfArrival(end)-1, 3);
    count=1;
    while ~isDone(bridge)
        truePosition(count, :)=bridge();
        count=count+1;
    end
    
    x=[truePosition(2:end-1, 1)', x(add_id:end)];
    y=[truePosition(2:end-1, 2)', y(add_id:end)];
else
    disp('error: please input the bridge_type!')
end

val_length = val_length + length(truePosition(2:end-1, 1));

%% Desired reference
q2=acos((x.^2+y.^2-l1^2-l2^2)/(2*l1*l2));

symb = 1;
for ij = 2:length(x)-1
    q2(ij) = symb * q2(ij);
    if q2(ij) == pi || q2(ij) == -pi ||  q2(ij) == 0 
        if sign(x(ij+1) - x(ij)) == sign(x(ij) - x(ij-1))
            symb = -symb;
        end
    end
end

q1=atan(y./x)-atan(l2*sin(q2)./(l1+l2*cos(q2)));

x1 = l1*cos(q1(1));
if round(x1, 2) ~= round(l1*cos(value_q(1)), 2)
    symbol_change_record = 1;
    q2 = -q2;
    q1=atan(y./x)-atan(l2*sin(q2)./(l1+l2*cos(q2)));
    for ij = 1:length(q1)
        if x(ij) < 0 && y(ij) > 0
            q1(ij)=q1(ij)+pi;
        elseif x(ij) < 0 && y(ij) < 0
            q1(ij)=q1(ij)+pi;
        elseif x(ij) > 0 && y(ij) < 0
            q1(ij)=q1(ij)+2*pi;
        end
    end
    
    for ij = 2:length(q1)
        if q1(ij)-q1(ij-1)<pi+0.1 && q1(ij)-q1(ij-1)>pi-0.1
            q1(ij)=q1(ij)-pi;
        end
        if q1(ij)-q1(ij-1)<-pi+0.1 && q1(ij)-q1(ij-1)>-pi-0.1
            q1(ij)=q1(ij)+pi;
        end
        if q1(ij)-q1(ij-1)>pi
            q1(ij:end)=q1(ij:end)-2*pi;
        end
        if q1(ij)-q1(ij-1)<-pi
            q1(ij:end)=q1(ij:end)+2*pi;
        end
        if q2(ij)-q2(ij-1)>pi
            q2(ij:end)=q2(ij:end)-2*pi;
        end
        if q2(ij)-q2(ij-1)<-pi
            q2(ij:end)=q2(ij:end)+2*pi;
        end
    end
end

wave_qdt=[(q1(2:val_length+1)-q1(1:val_length))', (q2(2:val_length+1)-q2(1:val_length))'];
wave_qdt=wave_qdt./dt;

q_control(2:val_length+1, :)=[q1(1:val_length)', q2(1:val_length)'];

qdt_control(2:val_length+1, :)=wave_qdt(1:val_length,:);

q2dt_control(1:val_length, :)=...
    (qdt_control(2:val_length+1,:)-qdt_control(1:val_length,:))/dt;

for ii = 1:val_length
    H11=m1*lc1^2+I1+m2*(l1^2+lc2^2+2*l1*lc2*cos(q_control(ii,2)))+I2;
    H12=m2*l1*lc2*cos(q_control(ii,2))+m2*lc2^2+I2;
    H21=H12;
    H22=m2*lc2^2+I2;
    h=m2*l1*lc2*sin(q_control(ii,2));
    
    part_1=-h*qdt_control(ii,2)*qdt_control(ii,1)-h*(qdt_control(ii,1)+qdt_control(ii,2))*qdt_control(ii,2);
    part_2=h*qdt_control(ii,1)*qdt_control(ii,1);
    
    tau_control(ii,1)=H11*q2dt_control(ii,1)+H12*q2dt_control(ii,2)+part_1;
    tau_control(ii,2)=H21*q2dt_control(ii,1)+H22*q2dt_control(ii,2)+part_2;
end

% plot_figures_inside = 1;
% if we need to visualization
if plot_figures_inside == 1
    figure();
    hold on
    plot(tau_control(2:val_length, 1))
    plot(tau_control(2:val_length, 2))
    ylabel('tau')
end

if plot_figures_inside == 1
    figure();
    hold on
    plot(qdt_control(1:val_length-1, 1), 'r')
    plot(qdt_control(1:val_length-1, 2), 'b')
    xlabel('time step')
    ylabel('dq/dt(control)')
    legend('dq/dt(1)', 'dq/dt(2)')
end

x_control=l1*cos(q_control(:,1))+l2*cos(q_control(:,1)+q_control(:,2));
y_control=l1*sin(q_control(:,1))+l2*sin(q_control(:,1)+q_control(:,2));

if plot_figures_inside == 1
    figure();
    hold on
    plot(x(1:val_length), y(1:val_length), 'r', 'LineWidth', 1)
    plot(x_control(1:val_length), y_control(1:val_length), 'b--', 'LineWidth', 2)
    xlim([-1, 1])
    ylim([-1, 1])
    legend('real', 'desired')
end

start_step = train_length -1000;
movie_step = 100;
time_all = train_length + 10000;

line_prop='solid';
% plot animation of the tracking control
if plot_movie == 1
    func_plot_movie(start_step, movie_step, time_all, q_control_all(:, 1), q_control_all(:, 2), properties, line_prop)
end

control_infor.q_control=q_control;
control_infor.qdt_control=qdt_control;
control_infor.q2dt_control=q2dt_control;
control_infor.tau_control=tau_control;

time_infor.val_length=val_length;

end

