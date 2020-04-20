clc
clear

%% 构造真实和测量数据
DATA_LENGTH = 30;%数据长度

rng();
w = randn(2,DATA_LENGTH);%生成DATA_LENGTH个标准正态分布的伪随机数 坐标的测量噪声
n = 0.25.*randn(2,DATA_LENGTH);%生成DATA_LENGTH个均值为0方差为0.5的正态分布的伪随机数 速度的测量噪声

p_truth = zeros(2,DATA_LENGTH);%生成位置的真实值
p_measured = zeros(2,DATA_LENGTH);%生成位置的测量值

v = zeros(2,DATA_LENGTH);%x方向上的初始速度为1m/s
v(:, 1) = [1,0].';
g = [0, 10];%y方向上重力加速度g

v_measured = zeros(2,DATA_LENGTH);
delta_t = 1;%每一次采样的时间间隔为1s

%生成真实值和测量值
for i = 1:DATA_LENGTH
    if i == 1
        v_measured(:,1) = v(:, 1) + w(1);%生成速度的测量值
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
    p_est(:, i) = p_est(:, i-1) + delta_t .* v_measured(:, i-1) + 0.5*power(delta_t, 2).*g.';%状态向量的传播
    sigma_propagate = A*sigma_propagate*A.' + Qt;%协方差矩阵的更新
    K = sigma_propagate * C.' * inv(C * sigma_propagate * C.' + R);%计算矩阵系数K
    p_est(:, i) = p_est(:, i) + K * (p_measured(:, i) - C * p_est(:, i));%更新测量得到的位置
    sigma_propagate = sigma_propagate - K * C * sigma_propagate;%更新协方差矩阵
end
%% 画图
plot(p_truth(1,:), p_truth(2,:), 'r');
hold on;
plot(p_est(1,:), p_est(2,:), 'b');



