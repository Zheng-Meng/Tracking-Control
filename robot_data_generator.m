function [xy, q, qdt, q2dt, tau] = robot_data_generator(time_infor, noise_level, dt, properties)

% generate time series for training and validation

[m1, m2, l1, l2, lc1, lc2, I1, I2] = matsplit(properties);

section_len = time_infor.section_len;
time_length = time_infor.time_length;

% generate noise as the control signal, note that we should smooth the data
% to make the generated signal continuous.
noise_interval=noise_level; % 3*10^(-2)
pert_length=time_length*2;
pert=-noise_interval+2*noise_interval*rand(pert_length,2);

pert=pert(100:end, :);

BB = smoothdata(pert, 'gaussian', 50);
pert = BB;

q=zeros(time_length,2);
qdt=zeros(time_length,2);
q2dt=zeros(time_length,2);
tau=zeros(time_length,2);
tau(:,:)=pert(1:time_length,:);

rng('shuffle')

% To avoid the values too large in random walk, every 'section_len' step we
% will reset the states of the two link robot arm.
for t_i = 1:time_length-1
    if mod(t_i, section_len)==0
        q(t_i, 1) = 2 * pi * rand(1);
        q(t_i, 2) = 2 * pi * rand(1) - pi;
        qdt(t_i,:)=[0,0];
    end
    
    H11=m1*lc1^2+I1+m2*(l1^2+lc2^2+2*l1*lc2*cos(q(t_i,2)))+I2;
    H12=m2*l1*lc2*cos(q(t_i,2))+m2*lc2^2+I2;
    H21=H12;
    H22=m2*lc2^2+I2;
    h=m2*l1*lc2*sin(q(t_i,2));
    
    part_1=-h*qdt(t_i,2)*qdt(t_i,1)-h*(qdt(t_i,1)+qdt(t_i,2))*qdt(t_i,2);
    part_2=h*qdt(t_i,1)*qdt(t_i,1);
    denominator=H12*H21-H11*H22;
    
    q2dt(t_i,1)=-(-part_1*H22+H12*part_2-H12*tau(t_i,2)+H22*tau(t_i,1))/denominator;
    q2dt(t_i,2)=-(part_1*H21-H11*part_2+H11*tau(t_i,2)-H21*tau(t_i,1))/denominator;
    
    if mod(t_i, section_len)==0
        q2dt(t_i,:)=[0,0];
        tau(t_i,:)=[0,0];
    end
    
    q(t_i+1,:)=q(t_i,:)+qdt(t_i,:)*dt;
    qdt(t_i+1,:)=qdt(t_i,:)+q2dt(t_i,:)*dt;
end

x=l1*cos(q(:,1))+l2*cos(q(:,1)+q(:,2));
y=l1*sin(q(:,1))+l2*sin(q(:,1)+q(:,2));
xy=[x, y];

end

