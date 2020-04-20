clc
clear

%% ������ʵ�Ͳ�������
DATA_LENGTH = 30;%���ݳ���

rng();
w = randn(2,DATA_LENGTH);%����DATA_LENGTH����׼��̬�ֲ���α����� ����Ĳ�������
n = 0.25.*randn(2,DATA_LENGTH);%����DATA_LENGTH����ֵΪ0����Ϊ0.5����̬�ֲ���α����� �ٶȵĲ�������

p_truth = zeros(2,DATA_LENGTH);%����λ�õ���ʵֵ
p_measured = zeros(2,DATA_LENGTH);%����λ�õĲ���ֵ

v = zeros(2,DATA_LENGTH);%x�����ϵĳ�ʼ�ٶ�Ϊ1m/s
v(:, 1) = [1,0].';
g = [0, 10];%y�������������ٶ�g

v_measured = zeros(2,DATA_LENGTH);
delta_t = 1;%ÿһ�β�����ʱ����Ϊ1s

%������ʵֵ�Ͳ���ֵ
for i = 1:DATA_LENGTH
    if i == 1
        v_measured(:,1) = v(:, 1) + w(1);%�����ٶȵĲ���ֵ
        continue;
    end
    
    p_truth(:, i) = p_truth(:, i-1) + delta_t .* v(:,i-1)+ 0.5*power(delta_t, 2).* g.';
    v(:,i) = (v(:,i-1) + delta_t.*g.').';
    p_measured(:, i) = p_truth(:, i) + w(i);
    v_measured(:, i) = v(:,i) + n(i);
end

%% Kalman filter
p_est = zeros(2,DATA_LENGTH);
% p_est(:, 1) = [5,5].'
A = eye(2);
C = eye(2);
sigma_propagate = diag([1,1]);
sigma_measured = eye(2);
Qt = diag([0.25, 0.25]);
R = diag([1,1]);

for i = 2:DATA_LENGTH
    p_est(:, i) = p_est(:, i-1) + delta_t .* v_measured(:, i-1) + 0.5*power(delta_t, 2).*g.';%״̬�����Ĵ���
    sigma_propagate = A*sigma_propagate*A.' + Qt;%Э�������ĸ���
    K = sigma_propagate * C.' * inv(C * sigma_propagate * C.' + R);%�������ϵ��K
    p_est(:, i) = p_est(:, i) + K * (p_measured(:, i) - C * p_est(:, i));%���²����õ���λ��
    sigma_propagate = sigma_propagate - K * C * sigma_propagate;%����Э�������
end
%% ��ͼ
plot(p_truth(1,:), p_truth(2,:), 'r');
hold on;
plot(p_est(1,:), p_est(2,:), 'b');



