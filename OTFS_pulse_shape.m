close all
clear all
clc

%delay Doppler grid dimensions
%span of delay axis
T = 50e-6;
%span of Doppler axis
Delf = 20e3;
%number of cells along Doppler
M = 20;
%Number of cells across Delay
N = 20;
%delay axis
delay_axis = 0:T/N:T-T/N;
%Doppler axis
doppler_axis = -Delf/2:Delf/N:Delf/2-Delf/N;

%% %transmit signal

dd_xmit = zeros(N,M);
%ideal pulse transmit signal in dd domain
%dd_xmit(1,2) = 1;
%dd_xmit(1,4) = 1;
%delay axis location of symbol
tx_symbol_n = [0];
%Doppler axis location of symbol
tx_symbol_m = [1];
for i = 1:length(tx_symbol_n)
    n=N/2-tx_symbol_n(i);m=M/2-tx_symbol_m(i);
    [M_ax,N_ax] = meshgrid((-M/2+m:M/2+m-1),(-N/2+n:N/2+n-1));
    dd_xmit = dd_xmit + sinc(M_ax/(0.5*pi)).*sinc(N_ax/(0.5*pi));
end

%% ISFFT to convert to time-frequency grid
X = fftshift(fft(fftshift(ifft(dd_xmit).',1)).',2);
%% Heisenberg transform to convert to time domain transmit signal
s_mat = ifft(X,M,1)*sqrt(M); % Heisenberg transform
s = s_mat(:);

%% channel each scatterer has delay and doppler coefficient

r = zeros(size(s));
channel_delay = [10,12,13];
channel_doppler = [2,3,4];
for i=1:length(channel_delay)
    r = r + circshift(s.*exp(-1j*2*pi*channel_doppler(i)/(M*N).*(1:M*N)'),-channel_delay(i));
end
% thermal noise at receiver
%snr in Db
signal_power = sum(s.*conj(s))/(M*N);
snr = 30;
sigma = signal_power/1e3;
r = awgn(r,snr);

time_axis = (0:1/(M*N):1-1/(M*N))*T;
%correct the range, there are M bins but range needs to be corrected
frequency_axis = 1:M; 

%% Receiver
r_mat = reshape(r,N,M);
% Wigner transform to convert to time frequency grid
Y = fft(r_mat,M,1)/sqrt(M); % Wigner transform
% SFFT to transform to DD domain again
dd_rx = fftshift(ifft(fftshift(fft(Y).')).');


%% plots
figure;
title("Transmit side");
subplot(3,1,1);
imagesc(doppler_axis,delay_axis, abs(dd_xmit))
subplot(3,1,2);
imagesc(doppler_axis,delay_axis, real(X))
subplot(3,1,3);
plot3(time_axis,abs(s).*cos(angle(s)),abs(s).*sin(angle(s)));
grid on


figure;
title("Receive side");
subplot(3,1,1);
plot3(time_axis,abs(r).*cos(angle(r)),abs(r).*sin(angle(r)));
grid on
subplot(3,1,2);
imagesc(doppler_axis,delay_axis, real(Y))
subplot(3,1,3);
imagesc(doppler_axis,delay_axis, abs(dd_rx))

% b = xcorr2(sinc(M_ax/(0.5*pi)).*sinc(N_ax/(0.5*pi)),dd_rx);
% imagesc(abs(b));
PAPR = max(s.*conj(s))/(sum(s.*conj(s))/length(s))


sine = sin(2*pi*time_axis);
PAPR = max(sine.*conj(sine))/(sum(sine.*conj(sine))/length(sine))

