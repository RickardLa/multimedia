%% TASK 1
% Step 1.1 Generate a stationary speech signal file
clc, clf, clear, close all

Fs = 12e3;            % Sampling frequency
Ts = 3;               % Duration of recording
nChannels= 1;         % Number of channels

% Set flag to 1 if new recording is neaded, else keep old
% WARNING!! Setting flag to 1 will overwrite old recording
recFlag = 0;          % Flag to make new recording

if recFlag == 1
    recObj = audiorecorder(Fs, 16, nChannels);   % Create audiorecorder object
    disp('START')
    recordblocking(recObj, Ts);                 % Record for Ts seconds
    disp('STOP') 
    s = getaudiodata(recObj);                   % Save audio as a double 
    audiowrite('Vowel.wav', s, Fs)              % Save as .wav in path
end

s = audioread('Vowel.wav');                 % read the audio-file
t = 0:1/Fs:(length(s)-1)/Fs;                % Convert samples to time

% Plotting time-domain signal
figure
% subplot(2,1,1)
plot(t,s)
grid on
xlabel('Time [s]')
title('Single tone')

% subplot(2,1,2)
% plot(t,s)
% grid on
% xlim([0 2])
% xlabel('Time [s]')
% title('Single tone zoomed in')

% subplot(3,1,3)
% plot(abs(fft(y)))
% grid on
% xlabel('Frequency [Hz]')
% title('FFT of single tone')

%% Step 1.2 Estimate the LPC model parameters
%clc, clf, clear, close all

p = 12;             % Model order
Ts = 0.3;           % Block size [s]
startSample = 8e3;  % What sample the block will start at 

[s, Fs] = audioread('Vowel.wav');                 % Read file from step 1.1 
y_block = s(startSample:startSample + Ts*Fs-1);   % Extract block of Ts seconds


% LPC on segment. 'A' contains model parameters, 'E' is the variance of the
% prediction errors 
[a, E] = lpc(y_block,p);
% KOLLA VAD SOM SKA SKE MED GAIN, G, tror det är vårt E vi får

%% Step 1.3 Calculate the residual sequence e(n)
t = 0:1/Fs:(length(y_block)-1)/Fs;


e = filter(a,1,y_block);              % Filter speech signal with A(z) to get residual sequence

figure
plot(t,e)
hold on
plot(t,y_block)
grid on
xlabel('Time [s]')
legend('Residual sequence', 'Original speech signal')


%% Step 1.4 Re-synthesize the speech using the estimated parameters

y_resynth = filter(1,a,e); 

figure
subplot(2,1,1)
plot(t,y_resynth)
title('Resynthesized speech')
xlabel('Time [s]')
grid on
subplot(2,1,2)
plot(t,y_block)
title('Original speech')
grid on
xlabel('Time [s]')



% ------------------------------------------------------------------------
% -------------------- TASK 2 -------------------------------------------
% ------------------------------------------------------------------------



%% Step 2.1 Record Sentence for 15 seconds
clc, clf, clear, close all

Fs = 12e3;              % Sampling frequency
Ts = 15;                % Duration of recording
nChannels= 1;           % Number of channels

% Set flag to 1 if new recording is neaded, else keep old
% WARNING!! Setting flag to 1 will overwrite old recording
recFlag = 0;          
if recFlag == 1
    recObj = audiorecorder(Fs, 16, nChannels);   % Create audiorecorder object
    disp('START')
    recordblocking(recObj, Ts);                 % Record for Ts seconds
    disp('STOP') 
    s = getaudiodata(recObj);                   % Save audio as a double 
    audiowrite('MySentence.wav', s, Fs)         % Save as .wav in path
end

s = audioread('MySentence.wav');                % read the audio-file
t = 0:1/Fs:(length(s)-1)/Fs;                    % Convert samples to time

% Plotting time-domain signal 
figure
plot(t,s)
grid on
xlabel('Time [s]')
title('15 second sentence')

%% STEP 2.2 Block based speech analysis
clc

p = 12;                             % number of parameters
Ts = 15;                            % Length of recorded signal [s]
durBlock = 0.02;                    % Duration of each block [s]
totBlocks = Ts/durBlock;            % Total number of blocks
L = length(s)/totBlocks;            % number of samples in each block [samples]
a = zeros(p+1, totBlocks);          % Parameter matrix a,
E = zeros(totBlocks,1);             % Variance matrix

% Compute each column of a in the same manner as in 1.2
% Loop through each block and find the parameters and put into matrix a
                      
for i=1:totBlocks
    y_block = s((i-1)*L+1:i*L);        % Extract the j-th block from y       
    [a(:,i), E(i)] = lpc(y_block,p);   % find parameters for current block
end




%% STEP 2.3 Block-based estimation of residual sequence e_hat

e_hat=zeros(totBlocks*L,1);           % Define residual vector

for i=1:totBlocks
    y_block = s((i-1)*L+1:i*L);                             % Extract the i-th block from y       
    e_hat((i-1)*L+1:i*L) = filter(a(:,i),1,y_block);        % Filter speech signal with A(z) to get residual sequence
end

t = 0:1/Fs:(length(s)-1)/Fs;                    % Convert samples to time

% Plotting time-domain signal 
figure
plot(t,e_hat)
grid on
xlabel('Time [s]')
title('Residuals for sentence')

% Plots the residual and compared with the signal
% Zoom in on two blocks (40ms and save)
figure
plot(t,e_hat)
hold on
plot(t,s)
grid on
xlabel('Time [s]')
legend('Residual sequence', 'Original speech signal')

%% STEP 2.4 Block-based speech re-synthesis

s_hat = zeros(length(s), 1);

for i=1:totBlocks
    y_block = s((i-1)*L+1:i*L);        % Extract the j-th block from y       
    s_hat((i-1)*L+1:i*L) = filter(1,a(:,i),e_hat((i-1)*L+1:i*L));
end

figure
subplot(3,1,1)
plot(t,s)
grid on
title( 'Original speech signal')
xlabel('Time [s]')

subplot(3,1,2)
plot(t,s_hat)
grid on
xlabel('Time [s]')
title('Re-synthesis sequence')

subplot(3,1,3)
plot(t,e_hat)
grid on
xlabel('Time [s]')
title('Residual sequence')




% ------------------------------------------------------------------------
% -------------------- TASK 3 -------------------------------------------
% ------------------------------------------------------------------------

%% Step 3.1 Form modified residual sequence e_tilde
% Use the residual sequence e_hat formed in 2.3
clc, close all, 
e_tilde = zeros(length(e_hat),1);
K = L/4;

for i=1:totBlocks
    sigVals = max(abs(e_hat((i-1)*L+1:L*i)),K);    % Find most significant vals
    [~,loc]=ismember(abs(e_hat),sigVals);       % Check positions
    indx = find(loc);                           % Extract non-zeros vals
    e_tilde(indx) = e_hat(indx);                % Extrac vals from hat to tilde
end

    
    
%% Step 3.2
% Repeat 2.3 with e_tilde

s_tilde = zeros(length(s), 1);

for i=1:totBlocks
    y_block = s((i-1)*L+1:i*L);        % Extract the j-th block from y       
    s_tilde((i-1)*L+1:i*L) = filter(1,a(:,i),e_tilde((i-1)*L+1:i*L));
end

figure
subplot(2,1,1)
plot(t,s)
legend( 'Original speech signal')
subplot(2,1,2)
plot(t,s_hat)
grid on
xlabel('Time [s]')
legend('Re-synthesis sequence')

writeFlag = 0;
if writeFlag == 1
    audiowrite('s_tilde.wav', s_tilde, Fs)
end


%%
for i=1:totBlocks
    y_block = s((i-1)*L+1:i*L);        % Extract the j-th block from y 
    y_block = y_block.*hamming(length(y_block));    % Multiply with hamming
    y = [y y_block zeros(length(y_block),1)] 
end




