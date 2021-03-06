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
% title('Single tone')

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

p = 0;            % Model order
Ts = 1;           % Block size [s]
startSample = 8e3;  % What sample the block will start at 

[s, Fs] = audioread('Vowel.wav');                 % Read file from step 1.1 
y_block = s(startSample:startSample + Ts*Fs-1);   % Extract block of Ts seconds


% LPC on segment. 'A' contains model parameters, 'E' is the variance of the
% prediction errors 
[a, E] = lpc(y_block,p);
 


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

p = 0;                             % number of parameters
Ts = 15;                           % Length of recorded signal [s]
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
% title('Residuals for sentence')

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
    s_hat((i-1)*L+1:i*L) = filter(1,a(:,i),e_hat((i-1)*L+1:i*L));
end

figure
subplot(3,1,1)
plot(t,s)
grid on
title('Original speech')
xlabel('Time [s]')

subplot(3,1,2)
plot(t,s_hat)
grid on
xlabel('Time [s]')
title('Re-synthesized speech')

subplot(3,1,3)
plot(t,e_hat)
grid on
xlabel('Time [s]')
title('Residual sequence')


audiowrite('shat.wav', s_hat, Fs)



% ------------------------------------------------------------------------
% -------------------- TASK 3 -------------------------------------------
% ------------------------------------------------------------------------

%% Step 3.1 Form modified residual sequence e_tilde
% Use the residual sequence e_hat formed in 2.3
clc, close all, 
e_tilde = zeros(length(e_hat),1);
K = 32;

for i=1:totBlocks
    sigVals = maxk(abs(e_hat((i-1)*L+1:L*i)),K);   % Find most significant vals. maxk() finds the k largest elements in e_hat.
    [~,loc]=ismember(abs(e_hat),sigVals);          % Check positions
    indx = find(loc);                              % Extract non-zeros vals
    e_tilde(indx) = e_hat(indx);                   % Extrac vals from e_hat to e_tilde
end

    
    
%% Step 3.2
% Repeat 2.3 with e_tilde

s_tilde = zeros(length(s), 1);

for i=1:totBlocks
    y_block = s((i-1)*L+1:i*L);        % Extract the i-th block from y       
    s_tilde((i-1)*L+1:i*L) = filter(1,a(:,i),e_tilde((i-1)*L+1:i*L));
end

figure
subplot(3,1,1)
plot(t,s)
grid on
xlabel('Time [s]')
title('Original speech')


subplot(3,1,2)
plot(t,s_tilde)
grid on
xlabel('Time [s]')
title('Re-synthesized speech')

subplot(3,1,3)
plot(t,e_tilde)
grid on
xlabel('Time [s]')
title('Modified residual sequence')




writeFlag = 1;
if writeFlag == 1
    audiowrite('stilde.wav', s_tilde, Fs)
end


% ------------------------------------------------------------------------
% -------------------- TASK 4 -------------------------------------------
% ------------------------------------------------------------------------

%% Step 4.1
clc, clf, clear, close

[s, Fs] = audioread('testPitch.wav');    % Read the audio-file

Ts = 8;                                  % Length of recorded signal [s]
durBlock = 0.02;                          % Duration of each block [s]
totBlocks = Ts/durBlock;                  % Total number of blocks
L = length(s)/totBlocks;                  % Number of samples in each block [samples]

padSize = 3;

cep1 = zeros((padSize+1)*L , totBlocks);
cep2 = zeros((padSize+1)*L , totBlocks);

for i=1:totBlocks
    x = s((i-1) * L + 1 : i * L);            % Divide original speech signal into blocks of length L
    x = x .* hamming(length(x));             % Multiply with Hamming-window to reduce Gibbs effect
    y = [x; zeros(padSize*length(x),1)];             % Pad with zeros
    
    cep1(:,i) = log10(abs(fft(y)));           % Spectrum in freq.
    cep2(:,i) = abs(ifft(cep1(:,i)));       % Spectrum in time
end


figure;
c = .5;
startBlock = 50;                    % Set start block between 1 and 220
t = (0:1/Fs:((padSize+1)*L-1)/Fs);


for j=startBlock:startBlock+20
    c2 = cep2(:,j) + j*c; 
    
    plot(t,c2)                        % Time domain plot
    xlabel('Time [s]')
    grid on
    hold on
%     title('Cepstrum for testPitch.wav')
    

    
end



% ------------------------------------------------------------------------
% -------------------- TASK 5 -------------------------------------------
% ------------------------------------------------------------------------
%% Step 5.1
clc, clf, clear, close all


[s_tilde, ~] = audioread('stilde.wav');    % Read the audio-file
[s, Fs] = audioread('MySentence.wav');    % Read the audio-file

p = 12; 
Ts = 15;                                  % Length of recorded signal [s]
durBlock = 0.02;                          % Duration of each block [s]
totBlocks = Ts/durBlock;                  % Total number of blocks
L = length(s)/totBlocks;                  % Number of samples in each block [samples]


N = 200; 
a_orig = zeros(p+1,totBlocks);
a_tilde = zeros(p+1, totBlocks); 

A_orig = zeros(N,totBlocks);
A_tilde = zeros(N,totBlocks);


for i=1:totBlocks
    x_orig = s((i-1) * L + 1 : i * L);
    x_tilde = s_tilde((i-1) * L + 1 : i * L);
    
    [a_orig(:,i), E_orig(i)] = lpc(x_orig,p);
    [a_tilde(:,i), E_tilde(i)] = lpc(x_tilde,p);
    
    A_orig(:,i) = 1./abs(fftshift(fft(a_orig(:,i),N))); 
    A_tilde(:,i) = 1./abs(fftshift(fft(a_tilde(:,i),N))); 
    
    d(i) = sum(10*log10(abs(A_orig(:,i)-A_tilde(:,i)).^2))/N;
end

f = Fs/2;   
w0 = 2*pi*f/Fs; 
w = (0:w0/(N-1):w0)/pi; 
figure;


idx = randi(totBlocks,1);
plot(w,pow2db((A_tilde(:,idx)).^2))
hold on
grid on
plot(w,pow2db((A_orig(:,idx)).^2))
legend('Re-synthesized signal', 'Original signal')
xlabel('pi*radians')
ylabel('dB')

figure;
plot(1:totBlocks, d)
xlabel('Block number')
ylabel('dB')
grid on



















