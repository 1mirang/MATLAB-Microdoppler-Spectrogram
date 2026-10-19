clc;
clearvars;
close all;
rng(1); % fix


%% 2. The 3D Geometrical Model
%Kích thước của căn phòng: X[-3, 3]m, Y[-3, 3]m, Z[0, 3]m

%Tọa độ phát và thu
Tx = [1.5, -1, 1.18];
Rx = [1.5, 1, 1.18];

%Tọa độ scatters
M = 6; %Số scatterers

SCF_x = -3 + 6 * rand(1,M);
SCF_y = -3 + 6 * rand(1,M);
SCF_z = 0 + 3 * rand(1,M);

%Con lắc
%Thông số của con lắc
L = 1.17;
xmax = 0.55;
g = 9.8;
dt = 0.001; % fix
t = 0:dt:15; %t=15s

%Quỹ đạo của con lắc theo 3 phương trình (38),(39),(40)
x_M = L * sin(asin(xmax/L) * cos(sqrt(g/L) * t));
y_M = zeros(1, length(t));
z_M = L * (1 - cos(asin(x_M/L)));

%%  VẼ ĐỒ THỊ MÔ PHỎNG MÔI TRƯỜNG 3D
figure('Name', 'Mô hình Không gian 3D - Micro-Doppler', 'Color', 'w');
hold on; grid on; view(3);
xlabel('Trục X (m)');
ylabel('Trục Y (m)');
zlabel('Trục Z (m)');
title('The 3D Geometrical Model');
axis([-3 3 -3 3 0 3]);

% Đầu thu phát
scatter3(Tx(1), Tx(2), Tx(3), 150, '^', 'DisplayName', 'Tx (máy phát)', 'MarkerFaceColor', 'r'); % label fix
scatter3(Rx(1), Rx(2), Rx(3), 150, 'v', 'DisplayName', 'Rx (máy thu)', 'MarkerFaceColor', 'b');  % label fix
% Fixed Scatterers
scatter3(SCF_x, SCF_y, SCF_z, 80, 's', 'DisplayName', 'Scatterer', 'MarkerFaceColor', [0.5 0.5 0.5] );
% Quỹ đạo Pendulum
plot3(x_M, y_M, z_M, 'g', 'LineWidth', 2.5, 'DisplayName', 'Quỹ đạo Pendulum');
% Starter point
scatter3(x_M(1), y_M(1), z_M(1), 100, 'o', 'MarkerFaceColor', 'g', 'DisplayName', 'Vị trí ban đầu của Pendulum');
legend('Location','best');


%% 3.Channel Transfer Function

% Tính khoảng cách dM_Tx (giữa Tx và con lắc) (10)
dM_Tx = sqrt((x_M - Tx(1)).^2 + (y_M - Tx(2)).^2 + (z_M - Tx(3)).^2);

% Tính khoảng cách dM_Rx (giữa con lắc và Rx) (11)
dM_Rx = sqrt((x_M - Rx(1)).^2 + (y_M - Rx(2)).^2 + (z_M - Rx(3)).^2);

%TV Propagation delay (16)
c0 = 299792458; % c0
tau_M = (dM_Tx + dM_Rx) / c0;


f0 = 5.32e9; %Carrier Freq
eta = 0.8; %Eta n
subc_f = 312.5e3; %Subcarrier frequencies difference

%OFDM subcarrier 
q = [-28:2:-2, -1, 1:2:27, 28]; % 30 sub
deltaf_q = q * subc_f; %(20)


%Chọn 6 scatterers cố định
M_fixed = M;

%Path gain (41)
cM_fixed = sqrt(2*eta);
cF_m = sqrt(2 * (1 - eta) / M_fixed);

%Random Phase Variable
theta_F = 2 * pi * rand(1, M_fixed) - pi; 
theta_M = 2 * pi * rand() - pi; 

%Tính CTF
H_F_m = sum(cF_m * exp(1i*theta_F)); %sigma 6 Scatterers %(19) % mean -> sum

H_total = zeros(length(q), length(t));
for idx = 1:length(q)
    H_M = cM_fixed * exp(1i * (theta_M - 2 * pi * (f0 + deltaf_q(idx)) * tau_M)); %(18)
    H_total(idx, :) = H_M + H_F_m; %(17)
end

%% 4.Spectrogram
% Tính tổng tín hiệu của 30 subcarriers (STFT tuyến tính nên tương đương tổng các STFT, Eq. 42)
H_signal = sum(H_total, 1); 
fs = 1/dt;               % Tần số lấy mẫu (1000 Hz do dt = 0.001)
win_length = 256;         % Số mẫu của window, lũy thừa của 2


%Gaussian window với sigma_w = 31.1 ms
sigma_t = 31.1e-3;
sigmasample = sigma_t / dt;                  % sigma tính theo số mẫu
alpha = (win_length - 1) / (2 * sigmasample); % gausswin: sigma = (N-1)/(2*alpha)
window = gausswin(win_length, alpha); 
noverlap = win_length - 5; % hop = 5 mẫu (5 ms)
nfft = 1024;             % Số điểm FFT

%Vẽ đồ thị
figure('Name', 'Microdoppler Spectrogram', 'Color', 'w');

% Hàm spectrogram
spectrogram(H_signal, window, noverlap, nfft, fs, 'yaxis', 'centered');

%Đồ thị
title('Spectrogram Microdoppler');
xlabel('Time, t (s)');
ylabel('Frequency, f (Hz)');
ylim([-50 50]);
%caxis([-15 25]);
h = gca;
mx = max(h.Children(1).CData(:));
clim([mx-50 mx])
colormap('hot');
