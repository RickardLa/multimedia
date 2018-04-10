%% Step 1.1 Generate a stationary speech signal file
clc
clf
clear 
close all

Fs = 12e3;            % Sampling frequency
Ts = 3;               % Duration of recording
nChannels= 1;         % Number of channels



recObj = audiorecorder(Fs, 8, nChannels);   % Create audiorecorder object
disp('START')
recordblocking(recObj, Ts);                 % Record for Ts seconds
disp('STOP') 

y = getaudiodata(recObj);                   % Save audio as a double 
audiowrite('Vowel.wav', y, Fs)              % Save as .wav in path
t = 0:1/Fs:(length(y)-1)/Fs;                % Convert samples to time



% Plotting time-domain signal 
subplot(2,1,1)
plot(t,y)
grid on
xlabel('Time [s]')
title('Single tone')

subplot(2,1,2)
plot(t,y)
grid on
xlim([0 2])
xlabel('Time [s]')
title('Single tone zoomed in')


% subplot(3,1,3)
% plot(abs(fft(y)))
% grid on
% xlabel('Frequency [Hz]')
% title('FFT of single tone')


%% Step 1.2 Estimate the LPC model parameters
clc
clf
clear 
close all

p = 12;             % Model order
Ts = 0.3;           % Block size [s]
startSample = 100;    % What sample the block will start at 

[y, Fs] = audioread('Vowel.wav');                 % Read file from step 1.1 
y_block = y(startSample:startSample + Ts*Fs-1);   % Extract block of Ts seconds


% LPC on segment. 'A' contains model parameters, 'E' is the variance of the
% prediction errors 
[a, E] = lpc(y_block,p)


%% Step 1.3 Calculate the residual sequence e(n)


t = 0:1/Fs:(length(y_block)-1)/Fs;
e = filter(1,a,y_block);              % Filter speech signal with A(z) to get residual sequence



plot(t,e)
hold on
plot(t,y_block)
grid on
xlabel('Time [s]')
legend('Residual sequence', 'Original speech signal')


%% Step 1.4 Re-synthesize the speech using the estimated parameters

y_resynth = filter(1,1/a(:),e);

plot(t,y_resynth)
hold on
plot(t,y_block)
grid on
xlabel('Time [s]')
legend('Re-synthesized speech', 'Original speech signal')

