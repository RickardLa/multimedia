
% %% Block 1 - Generate a Synthetic signal
% 
% clc, clear all, close all
% 
% t=[0:0.1:10*pi];
% signal = sin(t);
% 
% % figure
% % plot(signal)

%% Block 1
clc, clear, close all
originalImg = mat2gray(imread('lena.bmp','bmp')); % Load image file
[height, width] = size(originalImg); 

figure
imshow(originalImg)         % Display original image

R_comp = 0.5; % Compression ratio => Remove 256*256*0.5 coefficients
N1 = height-round(R_comp*16^2);  

% Step 2.2
L = 16;                                       % 16x16 pixels per block
totHeight = height/L;                         % Total number of height blocks
totWidth = width/L;                           % Total number of width blocks

vectorHeight = L * ones(1, totHeight);
vectorWidth = L * ones(1, totWidth);

% Create the blocks by converting image from a matrix to cell
allBlocks = mat2cell(originalImg, vectorHeight ,vectorWidth );

% Now compute DCT of all blocks and store them in DCTBlocks
DCTBlocks  = zeros(height, width);
compressedDCTBlocks  = zeros(height, width);
block = zeros(L); % Store the temporary block in.
compressedDCTBlock = zeros(L);  % Store temp block in for compression
TxVec = zeros(1,(height^2)*R_comp);    % Store all 1d koefficients here
n=1;    % Itterator for totVec

for i=1:totHeight       
    for j=1:totWidth   
        block = dct2(allBlocks{i,j});
        DCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L) = block; 
        
        % Find threshold value for each block (behövs ej?)
%         vectorBlock=reshape(block, 1, []); 
%         ascendDCT = sort(abs(vectorBlock),'descend');
%         thresholdBlock = ascendDCT(N1);
        
       % Performe zig-zag scan (set the last 128 values to zero)
       % Komprimmerar samtidigt som zig-zags and colt-45.
%        block = fliplr(triu(fliplr(block))); % Flipl and remove upper righ, flip back
%        for k=16:-1:9
%            block(k,16-(k-1))=0;
%        end
       
       blockVec = zigzag(block);    % Stolen algortithm , implement our own!!
       TxVec((n-1)*N1+1:N1*n) = blockVec(1:N1);        % compressed
       %totVecCheck((n-1)*width+1:width*n) = blockVec;   % Uncomprresed 
       
       %compressedDCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L)=block; 
       n = n+1; 
    end
end


%% Block 2 - Scalar quantification
% osaker pa detta steg, inte riktigt saker pa hur vardena och
% intervallen ska vara.. Det ska vara 8 bits => 256 nivåer, misstänker att
% längden på codebook ska vara 256 

% For Synthtetic signal with values in range from -1 to 1
% partition = linspace(-5,15,257); % To represent 257 intervalls, (will not take first and last)
% codebook = linspace(-5,15,256);
% 
% 
% [index,quants] = quantiz(TxVec,partition(2:end-1),codebook); % Quantize.
m = 8;      % Bits/sample

partition = linspace(min(TxVec),max(TxVec),2^m-1);
codebook = linspace(min(TxVec),max(TxVec),2^m); 
[index, quants] = quantiz(TxVec,partition, codebook);

% Plot signal and quantized signal
figure
plot(TxVec,'x')
hold on
plot(quants,'.')
legend('Original signal','Quantized signal');
axis([-.1 6.5 -1.2 1.2])
grid on



%% Block 3 Packetization

% PACKETIZATION
k = 127; % Symbols per packet
nPkts = ceil(length(index)/k); % Number of packets

%padd with zeros if last packet is not long enough
indexPad = zeros(1,nPkts*k);    
indexPad(1:length(index))=index;

% Create matrix with packages, number of rows equals number of packets
packetMatrix = reshape(indexPad,k,nPkts)';

%% Block 4 Reed Solomon encoding
m = 8; % Number of bits per symbol
n = 2^m-1;
% Block 9 Reed Solomon Decode
msgwords=gf(packetMatrix, m);
codes = rsenc(msgwords,n,k);
codewords = codes.x;

%% Block 9 Reed Solomon decoding
dec_msg = rsdec(codes,n,k);
dec_pktMtrx = dec_msg.x;

%% depacketization in Block-10
%dec_pktMtrx = packetMatrix;
[nPkts, nSyms] = size(dec_pktMtrx);
indexRx = reshape(dec_pktMtrx', 1,numel(dec_pktMtrx)); % Reshape to vector
indexRx = indexRx(1:256^2/2);

%% BLOCK 11 Generate the received DFT coefficients
%indexRx = index;
quantDCT = codebook(indexRx+1);

% Step 1.4 Verify
errorRatio = 1-(sum(quantDCT == quants))/length(quants) % Compute the error

% Plot stuff
figure
plot(TxVec)
hold on
plot(quantDCT)
axis([-.1 6.5 -1.2 1.2])
grid on
legend('Original signal','Received Quantized signal');

%% Block 12 Last block

% Put back zeros in order to reconstruct DCT matrix
newVec = zeros(1,width^2);

for i = 1:2:width*2
    j=ceil(i/2);
    newVec((i-1)*N1+1:N1*i) = quantDCT((j-1)*N1+1:N1*j);
end

RxDCTBlocks = cell(L);
colCount=0;
rowCount=1;

% PERFORM izigzag on all the sequences and put  in matrix
for i = 1:L^2:length(newVec)
    
    colCount=colCount+1;   
    if colCount==17
        rowCount=rowCount+1;
        colCount=1;
    end
    
    iBlock = izigzag(newVec(i:i+(L^2)-1),16,16);
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