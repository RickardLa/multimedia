% Block 1
clc, clear, close all
originalImg = mat2gray(imread('lena.bmp','bmp')); % Load image file
[height, width] = size(originalImg); 
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
totVec = zeros(1,(height^2)*R_comp);    % Store all 1d koefficients here
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
       block = fliplr(triu(fliplr(block))); % Flipl and remove upper righ, flip back
       for k=16:-1:9
           block(k,16-(k-1))=0;
       end
       blockVec = zigzag(block);    % Stolen algortithm , implement our own!!
       totVec((n-1)*N1+1:N1*n) = blockVec(1:N1);
       totVecCheck((n-1)*width+1:width*n) = blockVec;
       
       compressedDCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L)=block; 
       n = n+1;
    end
end



%% Last block
% input totVec, Output Lena

% Put back zeros in order to reconstruct DCT matrix

newVec = zeros(1,width^2);

for i = 1:2:width*2
    j=ceil(i/2);
    newVec((i-1)*N1+1:N1*i) = totVec((j-1)*N1+1:N1*j);
end

for i = 1:256
    for 
        newVec


