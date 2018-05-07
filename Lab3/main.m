
%% Block 1 - Generate a Synthetic signal

clc, clear all, close all

t=[0:0.1:10*pi];
signal = sin(t);

% figure
% plot(signal)

%% Block 1 - Generate signal EJ KLART; ANVÄND SYNTHETIC SIGNAL
% EJ FÄRDIGT!
% clc, clear, close all
% originalImg = mat2gray(imread('lena.bmp','bmp')); % Load image file
% [height, width] = size(originalImg); 
% imshow(originalImg)         % Display original image
% 
% R_comp = 0.5; % Compression ratio => Remove 256*256*0.5 coefficients
% N1 = height-round(R_comp*16^2);  
% 
% % Step 2.2
% L = 16;                                       % 16x16 pixels per block
% totHeight = height/L;                         % Total number of height blocks
% totWidth = width/L;                           % Total number of width blocks
% 
% vectorHeight = L * ones(1, totHeight);
% vectorWidth = L * ones(1, totWidth);
% 
% % Create the blocks by converting image from a matrix to cell
% allBlocks = mat2cell(originalImg, vectorHeight ,vectorWidth );
% 
% % Now compute DCT of all blocks and store them in DCTBlocks
% DCTBlocks  = zeros(height, width);
% compressedDCTBlocks  = zeros(height, width);
% compressedDCTBlock = zeros(1,L^2);  % Store temp block in for compression
% 
% for i=1:totHeight       
%     for j=1:totWidth   
%         block = dct2(allBlocks{i,j});
%         DCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L) = block; 
%         
%         % Find threshold value for each block
%         vectorBlock=reshape(block, 1, []); 
%         ascendDCT = sort(abs(vectorBlock),'descend');
%         thresholdBlock = ascendDCT(N1);
%         
%         % zig zag scan
%         
%         
%         % Compress
%         compressedDCTBlock = block;
%         compressedDCTBlock(abs(compressedDCTBlock)<=thresholdBlock) = 0;      
%         compressedDCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L) = compressedDCTBlock; 
%     end
% end
% 
% % Zig-Zag scanning - Generate 1d vector of DCT-coeffs
% 
% 
% %% Block 12 EJ FÄRDIGT
% % 2.4 Inverse DCT of compressedDCT
% % Create the blocks by converting image from a matrix to cell
% allCmprsdBlocks = mat2cell(compressedDCTBlocks, vectorHeight ,vectorWidth );
% 
% % Now compute IDCT of all blocks and store them in IDCTBlocks
% IDCTBlocks = zeros(height, width);
% for i=1:totHeight       
%     for j=1:totWidth   
%         block = idct2(allCmprsdBlocks{i,j});
%         IDCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L) = block; 
%     end
% end
% 
% %imwrite(IDCTBlocks, 'blockCompLena.bmp', 'bmp');
% errorImg = abs(originalImg - IDCTBlocks);
% errorImg30 = 30*abs(originalImg - IDCTBlocks);
% 
% PSNR = psnr(IDCTBlocks, originalImg)         % PSNR in dB
% SSIM = ssim(IDCTBlocks, originalImg)         % SSIM = 1 means the images are identical
% 
% colormap gray;
% subplot(2,2,1)
% imagesc(originalImg)
% title('Original')
% axis off;
% 
% subplot(2,2,2)
% imagesc(IDCTBlocks)
% title('Compressed')
% axis off;
% 
% subplot(2,2,3)
% imagesc(errorImg30)
% title('Error')
% axis off;
% 
% subplot(2,2,4)
% imshow(DCTBlocks), %colormap(gca,jet), colorbar
% title('DCT Domain')
% axis off;

%% Block 2 - Scalar quantification
% Valdigt osaker pa detta steg, inte riktigt saker pa hur vardena och
% intervallen ska vara.. Det ska vara 8 bits => 256 nivåer, misstänker att
% längden på codebook ska vara 256 , inte 8...

partition = [-1:.3:1]; % Length 7, to represent 7 intervals
codebook = [-1.2:.3:1]; % Length 8, one entry for each interval
[index,quants] = quantiz(signal,partition,codebook); % Quantize.

% Plot signal and quantized signal
plot(t,signal,'x',t,quants,'.')
legend('Original signal','Quantized signal');
axis([-.1 6.5 -1.2 1.2])
grid on



%% Block 3 Packetization

% PACKETIZATION
k = 127; % Symbols per packet
nPkts = ceil(length(index)/k); % Number of packets

% padd with zeros if last packet is not long enough
indexPad = zeros(1,nPkts*k);    
indexPad(1:length(index))=index;

% Create matrix with packages, number of rows equals number of packets
packetMatrix = reshape(indexPad,nPkts,k);

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
[nPkts, nSyms] = size(dec_pktMtrx);
indexRx = reshape(dec_pktMtrx, 1,numel(dec_pktMtrx)); % Reshape to vector

% Remove zero paddings
endIdx = find(indexRx == 0);    % find first padd
indexRx=indexRx(1:endIdx-1);      % Remove the zeros


%% BLOCK 11 Generate the received DFT coefficients
quantDCT = codebook(indexRx+1 );

% Step 1.4 Verify
errorRatio = 1-(sum(quantDCT == quants))/length(quants) % Compute the error

% Plot stuff
figure
plot(t,signal)
hold on
plot(t,quantDCT)
axis([-.1 6.5 -1.2 1.2])
grid on
legend('Original signal','Received Quantized signal');