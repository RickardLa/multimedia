%% TASK 1
clc, clear all, close all
% Step 1.1

t=0:0.1:10*pi;
signal = sin(t);

m = 8;      % Bits/sample

partition = linspace(min(signal),max(signal),2^m-1);
codebook = linspace(min(signal),max(signal),2^m); 
[idx, quant] = quantiz(signal,partition, codebook); 


quantRX = codebook(idx+1);

plot(signal)
hold on
plot(quantRX,'*')




%%
clc, clear all, close all

n = 10;
nc = 6;

A = randi(6,1,n)

B = reshape(A,[2 nc])
% 
% B = [B; zeros(nc,n/2)];
% 
% C = reshape(B, [1 n*2])