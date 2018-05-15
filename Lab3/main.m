%% Parameters
clc, clear, close all

R_comp = 0.5;                                 % Compression ratio
L = 16;                                       % 16x16 pixels per block
m = 8;                                        % Bits/sample
k = 127;                                      % Symbols per packet
errorRate = 0.2;                              % Lost packets [%]

%% Block 1
originalImg = mat2gray(imread('lena.bmp','bmp'));     % Load image file
[height, width] = size(originalImg);                  % Save dimensions

N1 = round(R_comp*L^2);                               % Number of coefficients removed
Nc = L^2-N1;                                          % Number of coefficients left

totHeight = height/L;                                 % Total number of height blocks
totWidth = width/L;                                   % Total number of width blocks

vectorHeight = L * ones(1, totHeight);
vectorWidth = L * ones(1, totWidth);

% Create the blocks by converting image from a matrix to cell
allBlocks = mat2cell(originalImg, vectorHeight ,vectorWidth );

% Now compute DCT of all blocks and store them in TxVec after zigzaging
block = zeros(L);                               % Store the temporary block
TxVec = zeros(1,floor((height^2)*R_comp));             % Store all coefficients in vector
n=1;                                            % Iterator for TxVec

for i=1:totHeight       
    for j=1:totWidth   
       block = dct2(allBlocks{i,j});
       blockVec = zigzag(block);                    % Stolen algortithm , implement our own!!
       TxVec((n-1)*Nc+1:Nc*n) = blockVec(1:Nc);     % Compressed
       n = n+1; 
    end
end


%% Block 2 - Scalar quantification
partition = linspace(min(TxVec),max(TxVec),2^m-1);
codebook = linspace(min(TxVec),max(TxVec),2^m); 
[index, quants] = quantiz(TxVec,partition, codebook);

% % Plot signal and quantized signal
% figure
% plot(TxVec,'x')
% hold on
% plot(quants,'.')
% legend('Original signal','Quantized signal');
% % axis([-.1 6.5 -1.2 1.2])
% grid on



%% Block 3 Packetization
nPkts = ceil(length(index)/k); % Number of packets

% Pad with zeros if last packet is not long enough
indexPad = zeros(1,nPkts*k);    
indexPad(1:length(index)) = index;

% Create matrix with packages, number of rows equals number of packets
packetMatrix = reshape(indexPad,k,nPkts)';

%% Block 4 Reed Solomon encoding
n = 2^m-1;
msgwords=gf(packetMatrix, m);
codes = rsenc(msgwords,n,k);


%% Block 5 Interleaving
clc
load('codewords.mat');
codewords = transpose(codes);
% codewords = codes;

%% Block 7 Bit errors and packet losses
clc

mode = "biterror";      % 'biterrors' or 'packetloss'
[row, col] = size(codewords);
switch mode
    case "biterror"
%         t = floor((n-k)/2);                                 % Number of erros per codeword
        t = 100;
        nw = length(codewords);                               % Number of codewords
        noise = randi(n,row,col) .* randerr(row,col,t);       % Generate noise-matrix
        noisyCode = codewords + noise;                        % Convert to same data type and add noise 
    case "packetloss"
        lostPackets = floor(0.2 * nPkts);               % Number of lost packets
        idx = randperm(n, lostPackets);                       % Non-repeating random indices
        codewords(idx,:) = 0;                                 % Set entire row to zero if lost
        noisyCode = codewords;
end


%% Block 8 De-interleaving
noisyCode = transpose(noisyCode); 
%% Block 9 Reed Solomon decoding
dec_msg = rsdec(noisyCode,n,k);
dec_pktMtrx = dec_msg.x;

%% Block-10 de-packetization
[nPkts, nSyms] = size(dec_pktMtrx);
indexRx = reshape(dec_pktMtrx', 1,numel(dec_pktMtrx)); % Reshape to vector
indexRx = indexRx(1:height^2*R_comp);



%% BLOCK 11 Generate the received DFT coefficients

quantDCT = codebook(indexRx+1);

% Step 1.4 Verify
errorRatio = 1-(sum(quantDCT == quants))/length(quants) % Compute the error

% Plot stuff
% figure
% plot(TxVec)
% hold on
% plot(quantDCT)
% axis([-.1 6.5 -1.2 1.2])
% grid on
% legend('Original signal','Received Quantized signal');

%% Block 12 Last block

% Put back zeros in order to reconstruct DCT matrix

tmp = reshape(quantDCT, [Nc width]);
tmp = [tmp; zeros(Nc, width)];
newVec = reshape(tmp, [1 width^2]);


RxDCTBlocks = cell(L);
colCount=0;
rowCount=1;

% PERFORM izigzag on all the sequences and put  in matrix
for i = 1:L^2:length(newVec)
    
    colCount=colCount+1;   
    if colCount == L+1
        rowCount=rowCount+1;
        colCount=1;
    end
    
    iBlock = izigzag(newVec(i:i+(L^2)-1),L,L);
    RxDCTBlocks{rowCount,colCount} = iBlock;
  
end    

% Now compute IDCT of all blocks and store them in IDCTBlocks
IDCTBlocks = zeros(height, width);
for i=1:totHeight       
    for j=1:totWidth
        block = idct2(RxDCTBlocks{i,j});
        IDCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L) = block; 
        
    end
end



figure
colormap gray
imshow(IDCTBlocks)
title('Compressed')
axis off;

PSNR = psnr(IDCTBlocks, originalImg)         % PSNR in dB
SSIM = ssim(IDCTBlocks, originalImg)         % SSIM = 1 means the images are identical