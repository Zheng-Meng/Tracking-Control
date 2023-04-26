function [Wout, r_end] = func_reservoir_train(data_reservoir, time_infor, input_infor, res_infor, dim_in, dim_out)
rng('shuffle');
% train the reservoir computing by given the input and output.
xy=data_reservoir.xy;
q=data_reservoir.q;
qdt=data_reservoir.qdt;
q2dt=data_reservoir.q2dt;
tau=data_reservoir.tau;

washup_length=time_infor.washup_length;
train_length=time_infor.train_length; % 50000

W_in=res_infor.W_in;
res_net=res_infor.res_net;
alpha=res_infor.alpha;
kb=res_infor.kb;
beta=res_infor.beta;
n=res_infor.n;

r_train=zeros(n,train_length-washup_length);
y_train=zeros(dim_out,train_length-washup_length);

train_x=zeros(train_length,dim_in);
train_y=zeros(train_length,dim_out);
if length(input_infor)==2 && strcmp(input_infor(1), 'xy') == 1 && strcmp(input_infor(2), 'qdt') == 1
    train_x(:,:)=[xy(1:train_length,:),xy(2:train_length+1,:), qdt(1:train_length,:),qdt(2:train_length+1,:)];
elseif length(input_infor)==1 && strcmp(input_infor(1), 'q') == 1
    train_x(:,:)=[q(1:train_length,:),q(2:train_length+1,:)];
elseif length(input_infor)==2 && strcmp(input_infor(1), 'q') == 1 && strcmp(input_infor(2), 'qdt') == 1
    train_x(:,:)=[q(1:train_length,:),q(2:train_length+1,:), qdt(1:train_length,:),qdt(2:train_length+1,:)];
elseif length(input_infor)==3 && strcmp(input_infor(1), 'xy') == 1 && strcmp(input_infor(2), 'qdt') == 1 && strcmp(input_infor(3), 'q2dt') == 1 
    train_x(:,:)=[q(1:train_length,:),q(2:train_length+1,:), qdt(1:train_length,:),qdt(2:train_length+1,:), q2dt(1:train_length,:),q2dt(2:train_length+1,:)];
end

train_y(:,:)=tau(1:train_length,:);
train_x=train_x';
train_y=train_y';

r_all=zeros(n,train_length+1);%2*rand(n,1)-1;%
for ti=1:train_length
    r_all(:,ti+1)=(1-alpha)*r_all(:,ti) + alpha*tanh(res_net*r_all(:,ti)+W_in*train_x(:,ti)+kb*ones(n,1));
end

r_out=r_all(:,washup_length+2:end); % n * (train_length - 11)
r_out(2:2:end,:)=r_out(2:2:end,:).^2;
r_end(:)=r_all(:,end); % n * 1

r_train(:,:) = r_out;
y_train(:,:) = train_y(1:dim_out,washup_length+1:end);
% linear regression
Wout=y_train*r_train'*(r_train*r_train'+beta*eye(n))^(-1);

end

