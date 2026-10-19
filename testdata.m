clc; 
clearvars; 
close all;

%NẠP DỮ LIỆU
load('RL_grid_corridor.mat'); 

% Chọn ra 1 điểm đo bất kì với 4310 điểm lấy mẫu
RL_column = size(RL, 2); %size(...,2) lấy số cột (dim2)
rand_idx = randi(RL_column);

% Lấy ngẫu nhiên 1 cột và lật ngang thành mảng 1x4310
H_signal = RL(:, rand_idx).'; 

% Chuẩn hóa biên độ Linfinity
H_signal = H_signal / max(abs(H_signal)); 

%% SPECTROGRAM SETUP
Ts = 2.6667e-11;          
fs = 1 / Ts;
win_length = 256;     
noverlap = win_length - 1;
nfft = 2048; %2048
% Window (hamming, hann, blackman, kaiser)
window = hamming(win_length); 

%% SPECTROGRAM
figure('Name', 'Fast-Time Spectrogram', 'Color', 'w');

% Vẽ spectrogram
spectrogram(H_signal, window, noverlap, nfft, fs, 'yaxis', 'centered');

%Đồ thị
title(sprintf('Realtime Data Spectrogram - Cột: %d', rand_idx));
xlabel('Time, t (s)');
ylabel('Frequency, f (Hz)');
%clim([max_dB - 40, max_dB]);
colormap('hot');
colorbar;
