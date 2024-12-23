% 清除工作区中的所有变量，关闭所有图形窗口，清除命令窗口中的所有内容，
% 并关闭所有警告信息的显示。
clear; close all; clc; warning off;

%% 1. 抽样模拟信号

% 定义模拟信号的连续时间采样率和时间向量
fs_cont = 0.01;  % 连续时间采样间隔 (s)
t_cont = 0:fs_cont:10;  % 时间轴，仿真信号时间长度为10秒
% 定义模拟信号：由3个正弦波叠加而成
signal_cont = 0.1*cos(0.15*pi*t_cont) + 1.5*sin(2.5*pi*t_cont) + 0.5*cos(4*pi*t_cont);

% 进行信号的抽样，定义抽样频率和时间向量
fs = 20;  % 抽样频率为20Hz
fs_sample = 1/fs;  % 抽样间隔（s）
t_sample = 0:fs_sample:10;  % 抽样信号的时间轴
% 生成抽样后的信号
signal_sample = 0.1*cos(0.15*pi*t_sample) + 1.5*sin(2.5*pi*t_sample) + 0.5*cos(4*pi*t_sample);

%% 2. 计算频谱

% 原始模拟信号的频谱计算
n_cont = length(signal_cont);  % 信号长度
f_cont = fs_cont * (-n_cont/2:n_cont/2-1) / n_cont;  % 频率轴，重新定义频率范围
signal_cont_f = fftshift(fft(signal_cont));  % 对信号进行快速傅里叶变换 (FFT) 并调整频谱的零频位置

% 抽样信号的频谱计算
n_sample = length(signal_sample);  % 抽样信号的长度
f_sample = fs * (-n_sample/2:n_sample/2-1) / n_sample;  % 频率轴，重新定义频率范围
signal_sample_f = fftshift(fft(signal_sample));  % 对抽样信号进行FFT并调整频谱

%% 3. 绘制时域和频域图形

% 原始模拟信号时频域图
hFig = figure;
set(hFig, 'Name', '1.原始信号时频域', 'NumberTitle', 'off');
subplot(2,1,1);
plot(t_cont, signal_cont);
title('原模拟信号');
xlabel('时间 (s)');
ylabel('信号幅度');

subplot(2,1,2);
plot(f_cont, abs(signal_cont_f));
title('原模拟信号频谱');
xlabel('频率 (Hz)');
ylabel('幅度');
exportgraphics(gcf, 'results/1.原始信号时频域.png', 'Resolution', 300);

% 抽样信号时频域图
hFig = figure;
set(hFig, 'Name', ['2.抽样信号时频域（' num2str(fs) 'Hz采样）'], 'NumberTitle', 'off');
subplot(2,1,1);
stem(t_sample, signal_sample, 'r');
title('抽样信号');
xlabel('时间 (s)');
ylabel('信号幅度');

subplot(2,1,2);
plot(f_sample, abs(signal_sample_f));
title('抽样信号频谱');
xlabel('频率 (Hz)');
ylabel('幅度');
exportgraphics(gcf, ['results/2.抽样信号时频域（' num2str(fs) 'Hz采样）.png'], 'Resolution', 300);

%% 4. A律编码与解码

% A律编码参数设置
A = 87.6;  % A律编码的参数

% 归一化处理，使信号的幅度在[-1, 1]之间
signal_sample = signal_sample / max(abs(signal_sample));

% A律编码：对抽样信号进行A律编码
encoded_signal = A_law_encoding(signal_sample, A);

% 转换为8位PCM码
encoded_8bit = round((encoded_signal + 1) * 127);

% A律解码：对编码后的信号进行解码，恢复原始信号
decoded_signal = A_law_decoding(encoded_signal, A);

% 计算量化误差
quantization_error = signal_sample - decoded_signal;

% 计算均方误差 (MSE)
MSE = mean(quantization_error.^2);

% 输出编码信息和量化误差
disp('8位编码：');
disp(encoded_8bit);
disp('译码值：');
disp(decoded_signal);
disp('量化误差：');
disp(quantization_error);
disp(['均方误差(MSE): ', num2str(MSE)]);

%% 5. 绘制编码信号与译码信号对比图

hFig = figure;
set(hFig, 'Name', '3.编码信号与译码信号对比图', 'NumberTitle', 'off');
subplot(2,1,1);
plot(t_sample, encoded_signal, 'b', 'LineWidth', 1.5);
title('A律编码后的信号');
xlabel('时间 (秒)');
ylabel('编码信号幅度');

subplot(2,1,2);
plot(t_sample, decoded_signal, 'r', 'LineWidth', 1.5);
title('译码信号');
xlabel('时间 (秒)');
ylabel('幅度');
exportgraphics(gcf, 'results/3.编码信号与译码信号.png', 'Resolution', 300);

%% 6. 绘制量化误差图

hFig = figure;
set(hFig, 'Name', '4.量化误差图', 'NumberTitle', 'off');
plot(t_sample, quantization_error, 'g', 'LineWidth', 1.5);
title('量化误差');
xlabel('时间 (秒)');
ylabel('误差幅度');
exportgraphics(gcf, 'results/4.量化误差图.png', 'Resolution', 300);

%% 7. 输出结果到 Excel 文件

% 创建包含编码、解码和量化误差数据的表格
T = table(t_sample', encoded_8bit', decoded_signal', quantization_error', ...
    'VariableNames', {'Time', '8位编码', '译码值', '量化误差'});

% 将表格保存到 Excel 文件
filename = 'results/编译码数据.xlsx';
writetable(T, filename);

disp(['数据已保存到 Excel 文件: ', filename]);

%% 过程中的函数：A 律编解码

% A 律编码函数
function encoded_signal = A_law_encoding(signal, A)
    % A 律编码公式
    encoded_signal = sign(signal) .* log(1 + A * abs(signal)) / log(1 + A);
end

% A 律解码函数
function decoded_signal = A_law_decoding(encoded_signal, A)
    % A 律解码公式
    decoded_signal = sign(encoded_signal) .* ((1 + A).^abs(encoded_signal) - 1) / A;
end
