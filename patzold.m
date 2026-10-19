%% PATZOLD_FIXED.M
% Mo phong Micro-Doppler cua con lac (Abdelgawwad, Borhani, Paetzold - Sensors 2020)
% Cac thay doi so voi ban goc duoc danh dau [SUA #n]
clc; clearvars; close all;
rng(1);                                   % [SUA #1] seed dat TRUOC moi lenh rand

%% 1. Tham so
c0  = 299792458;                          % [SUA #7] khong can toolbox (physconst)
f0  = 5.32e9;                             % tan so song mang
eta = 0.8;
M   = 6;                                  % so fixed scatterer
df  = 312.5e3;                            % khoang cach subcarrier
q   = [-28:2:-2, -1, 1:2:27, 28];         % [SUA #4] 30 subcarrier (Eq. 20, Intel 5300)

Tx = [1.5, -1, 1.18];
Rx = [1.5,  1, 1.18];

L = 1.17;  xmax = 0.55;  g = 9.8;         % con lac
fs = 1000;                                % [SUA #5] 1 kHz nhu CSI trong paper (f_D,max ~ 40 Hz)
dt = 1/fs;
t  = 0:dt:15;

%% 2. Hinh hoc 3D
SCF = [-3+6*rand(M,1), -3+6*rand(M,1), 3*rand(M,1)];   % chi dung de ve hinh

x_M = L*sin(asin(xmax/L)*cos(sqrt(g/L)*t));            % Eq. (38)
y_M = zeros(size(t));                                  % Eq. (39)
z_M = L*(1 - cos(asin(x_M/L)));                        % Eq. (40)

figure('Name','Mo hinh 3D','Color','w'); hold on; grid on; view(3);
scatter3(Tx(1),Tx(2),Tx(3),150,'^','filled','MarkerFaceColor','r','DisplayName','Tx (may phat)');
scatter3(Rx(1),Rx(2),Rx(3),150,'v','filled','MarkerFaceColor','b','DisplayName','Rx (may thu)');   % [SUA #6]
scatter3(SCF(:,1),SCF(:,2),SCF(:,3),80,'s','filled','MarkerFaceColor',[0.5 0.5 0.5],'DisplayName','Fixed scatterer');
plot3(x_M,y_M,z_M,'g','LineWidth',2.5,'DisplayName','Quy dao con lac');
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis([-3 3 -3 3 0 3]);
title('The 3D Geometrical Model');
legend('Location','best');                % [SUA #6] legend SAU khi ve

%% 3. Channel Transfer Function
dTx = sqrt((x_M-Tx(1)).^2 + (y_M-Tx(2)).^2 + (z_M-Tx(3)).^2);   % Eq. (10)
dRx = sqrt((x_M-Rx(1)).^2 + (y_M-Rx(2)).^2 + (z_M-Rx(3)).^2);   % Eq. (11)
tau = (dTx + dRx)/c0;                                           % Eq. (16)

cM  = sqrt(2*eta);                        % Eq. (41)
cF  = sqrt(2*(1-eta)/M);
thM = 2*pi*rand - pi;
thF = 2*pi*rand(1,M) - pi;

HF = sum(cF*exp(1i*thF));                 % [SUA #3] TONG theo Eq. (17), khong phai mean

fq  = f0 + q(:)*df;                                   % 30 x 1
HM  = cM*exp(1i*(thM - 2*pi*fq*tau));                 % Eq. (18), kich thuoc 30 x numel(t)
H   = HM + HF;                                        % Eq. (17)
Hs  = sum(H,1);                                       % tong STFT theo q = STFT cua tong (tuyen tinh), Eq. (42)

%% 4. Spectrogram
sigma_t = 31.1e-3;                        % sigma_w trong paper
Lw      = 256;
sigma_n = sigma_t*fs;                     % do lech chuan theo so mau
alpha   = Lw/(2*sigma_n);                 % [SUA #2] gausswin: sigma = Lw/(2*alpha) mau
win     = gausswin(Lw, alpha);
hop     = 5;                              % buoc truot (giam khoi luong tinh toan so voi hop = 1)
noverlap= Lw - hop;
nfft    = 1024;

[S, f, tt] = spectrogram(Hs, win, noverlap, nfft, fs, 'centered');
P    = abs(S).^2;
PdB  = 10*log10(P/max(P(:)));             % [SUA #8] dB tuong doi so voi dinh

figure('Name','Micro-Doppler Spectrogram','Color','w');
imagesc(tt, f, PdB); axis xy;
xlabel('Time, t (s)'); ylabel('Frequency, f (Hz)');
title('Spectrogram Micro-Doppler (mo phong)');
ylim([-60 60]);
clim([-40 0]);                            % R2022a tro len; ban cu dung caxis
colormap('hot'); colorbar;

%% 5. Doi chung dinh luong: mean Doppler shift (Eq. 25 va 37)
fM    = -f0*gradient(tau, dt);            % Doppler cua con lac, f_M = -f0 * d(tau)/dt
B_ana = cM^2*fM/(cM^2 + M*cF^2);          % Eq. (25)
B_spec = (f(:).'*P)./sum(P,1);            % Eq. (37)
B_ref  = interp1(t, B_ana, tt(:).', 'linear');
nmse   = sum((B_spec - B_ref).^2)/sum(B_ref.^2);

figure('Name','Mean Doppler','Color','w'); hold on; grid on;
plot(t, B_ana, 'b', 'LineWidth', 1.5, 'DisplayName', 'B giai tich (Eq. 25)');
plot(tt, B_spec, 'r--', 'LineWidth', 1.2, 'DisplayName', 'B tu spectrogram (Eq. 37)');
xlabel('Time, t (s)'); ylabel('Mean Doppler shift (Hz)');
title(sprintf('NMSE = %.4f', nmse)); legend('Location','best');
