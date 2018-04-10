%% Step 1.1
clc
clf
clear 
close all

Fs = 12e3;            % Sampling frequency
Ts = 3;               % Duration of recording
nChannels= 1;         % Number of channels

% hej

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


%% Step 1.2
clc
clf
clear 
close all

p = 12;             % Model order
Ts = 0.3;           % Block size [s]
startSample = 1;    % What sample the block will start at 

[y, Fs] = audioread('Vowel.wav');           % Read file from step 1.1 
y = y(startSample:startSample + Ts*Fs-1);   % Extract block of Ts seconds
t = 0:1/Fs:(length(y)-1)/Fs;                % Convert samples to time


% plot(t,y)
% xlabel('Time [s]')

% LPC on segment. 'A' contains model parameters, 'E' is the variance of the
% prediction errors 
[A, E] = lpc(y,p)


%% Step 1.3
clc
clf
clear 
close all



