% ------------------------------------------------------------------------
% -------------------- TASK 1.1 ------------------------------------------
% ------------------------------------------------------------------------

%% Step 1.1
clc, clear, close

video = VideoReader('Trees1.avi');
width = video.Width;
height = video.Height;


mov = struct('frames',zeros(height,width));    % Struct to save frames in

numberOfFrames = round(video.Duration * video.FrameRate);

for i=1:numberOfFrames
    mov(i).frames = readFrame(video);
end

% Saving two Grayscale-images, 10 frames apart. 
I1 = rgb2gray(mov(10).frames);
I11 = rgb2gray(mov(11).frames); 
I2 = rgb2gray(mov(20).frames);

% imwrite(I1, 'Figures/frame10.bmp', 'bmp');
% imwrite(I11, 'frame11.bmp', 'bmp');
% imwrite(I2, 'Figures/frame20.bmp', 'bmp');

imagesc(I1), colormap gray
figure
imagesc(I2), colormap gray

%% Step 1.2
clc, close

originalImg = mat2gray(I1); 
DCTCoeff = dct2(originalImg); 


imshow(((DCTCoeff))), %colormap(gca,jet), colorbar


%% Step 1.3
clc, close

th = floor(0.9 * width * height);

% Reshape matrix to 1-column and order DCTCoeff in ascending order
vectorDCT = reshape(DCTCoeff, 1, []); 
ascendDCT = sort(abs(vectorDCT)); 

% Assign threshold 
threshold = ascendDCT(th);

% If abs(coefficient)<= threshold, set to 0
compressedDCT = DCTCoeff;
compressedDCT(abs(compressedDCT)<=threshold) = 0; 

%% Step 1.4
clc, close

% Inverse DCT of compressedDCT
compressedImg = idct2(compressedDCT);

errorImg = abs(originalImg - compressedImg);

PSNR = psnr(compressedImg, originalImg)         % PSNR in dB
SSIM = ssim(compressedImg, originalImg)         % SSIM = 1 means the images are identical

colormap gray;
subplot(2,2,1)
imagesc(originalImg)
title('Original')
axis off;

subplot(2,2,2)
imagesc(compressedImg)
title('Compressed')
axis off;

subplot(2,2,3)
imagesc(errorImg)
title('Error')
axis off;

subplot(2,2,4)
imshow(((DCTCoeff))), %colormap(gca,jet)
axis off;

%saveas(gca, 'Figures/Compressed.eps','epsc');

%% ------------------------------------------------------------------------
% -------------------- TASK 1.2 -------------------------------------------
% ------------------------------------------------------------------------
clc,clear,close;

% Step 2.1
originalImg = mat2gray(imread('Figures/frame10.bmp','bmp'));
[height, width] = size(originalImg); 

% Removing border-pixels so (height*width)/blockSize has no remainder
height = height - 1; 
width = width - 1; 
originalImg = originalImg(1:height,1:width);     

% Step 2.2
L = 8;                                        % 8x8 pixels per block
totHeight = height/L;                         % Total number of height blocks
totWidth = width/L;                           % Total number of width blocks

vectorHeight = L * ones(1, totHeight);
vectorWidth = L * ones(1, totWidth);

% Create the blocks by converting image from a matrix to cell
allBlocks = mat2cell(originalImg, vectorHeight ,vectorWidth );

% Now compute DCT of all blocks and store them in DCTBlocks
DCTBlocks = zeros(height, width);
for i=1:totHeight       
    for j=1:totWidth   
        DCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L) = dct2(allBlocks{i,j}); 
    end
end

% 2.3 COMPRESS
th = floor(0.9 * width * height);

% Reshape matrix to 1-column and order DCTCoeff in ascending order
vectorDCTBlocks = reshape(DCTBlocks, 1, []); 
ascendDCT = sort(abs(vectorDCTBlocks)); 

% Assign threshold 
thresholdBlocks = ascendDCT(th);

% If abs(coefficient)<= threshold, set to 0
compressedDCTBlocks = DCTBlocks;
compressedDCTBlocks(abs(compressedDCTBlocks)<=thresholdBlocks) = 0; 

% 2.4 Inverse DCT of compressedDCT
% Create the blocks by converting image from a matrix to cell
allCmprsdBlocks = mat2cell(compressedDCTBlocks, vectorHeight ,vectorWidth );

% Now compute IDCT of all blocks and store them in IDCTBlocks
IDCTBlocks = zeros(height, width);
for i=1:totHeight       
    for j=1:totWidth   
        block = idct2(allCmprsdBlocks{i,j});
        IDCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L) = block; 
    end
end

imwrite(IDCTBlocks, 'blockComp10.bmp', 'bmp');
errorImg = abs(originalImg - IDCTBlocks);
errorImg30 = 30*abs(originalImg - IDCTBlocks);

PSNR = psnr(IDCTBlocks, originalImg)         % PSNR in dB
SSIM = ssim(IDCTBlocks, originalImg)         % SSIM = 1 means the images are identical

colormap gray;
subplot(2,2,1)
imagesc(originalImg)
title('Original')
axis off;

subplot(2,2,2)
imagesc(IDCTBlocks)
title('Compressed')
axis off;

subplot(2,2,3)
imagesc(errorImg30)
title('Error')
axis off;

subplot(2,2,4)
imshow(DCTBlocks), %colormap(gca,jet), colorbar
title('DCT Domain')
axis off;


%% ------------------------------------------------------------------------
% -------------------- TASK 1.3 -------------------------------------------
% ------------------------------------------------------------------------
clc, 
waveletCompressedImg = mat2gray(imread('Figures/waveletCompressed10.bmp','bmp'));
originalImg = mat2gray(imread('Figures/frame10.bmp','bmp'));

errorImg = abs(originalImg - waveletCompressedImg);
errorImg30 = 30*abs(originalImg - waveletCompressedImg);

PSNR = psnr(waveletCompressedImg, originalImg)         % PSNR in dB
SSIM = ssim(waveletCompressedImg, originalImg)         % SSIM = 1 means the images are identical

figure
colormap gray;
subplot(2,2,1)
imagesc(originalImg)
title('Original')
axis off;

subplot(2,2,2)
imagesc(waveletCompressedImg)
title('Compressed Wave')
axis off;

subplot(2,2,3)
imagesc(errorImg30)
title('Error: 30*(Original - Wave)')
axis off;

subplot(2,2,4)
imagesc(IDCTBlocks), %colormap(gca,jet)
title('Compressed DCT Blocks')
axis off;


%% ------------------------------------------------------------------------
% -------------------- B TASK 2.2 -------------------------------------------
% ------------------------------------------------------------------------
clc, clear, close all
Inew = mat2gray(imread('frame11.bmp','bmp'));
%Iold = mat2gray(imread('blockComp10.bmp','bmp'));
Iold = mat2gray(imread('waveletCompressed10.bmp','bmp'));
Iold_orig = mat2gray(imread('Figures/frame10.bmp','bmp'));

% Removing border-pixels so (height*width)/blockSize has no remainder
Inew = Inew(1:end-1,1:end-1);
Iold = Iold(1:end-1,1:end-1);

th=45/255;
Idiff = abs(Iold-Inew);
Idiff(abs(Idiff) < th)=0;

[height, width] = size(Idiff); 
height = height -1; 
width = width -1;
% 
% Idiff = Idiff(1:height,1:width);
% Iold = Iold(1:height,1:width);
% Inew = Inew(1:height,1:width);
L = 16;                                       % 16x16 pixels per block
totHeight = ceil(height/L);                         % Total number of height blocks
totWidth = ceil(width/L);                           % Total number of width blocks

vectorHeight = L * ones(1, totHeight);
vectorWidth = L * ones(1, totWidth);

IdiffCell = mat2cell(Idiff, vectorHeight, vectorWidth);
oldCell = mat2cell(Iold, vectorHeight, vectorWidth);
newCell = mat2cell(Inew, vectorHeight, vectorWidth);

%% MOTION MATRIX

Imotion = zeros(height,width);

for i = 1:totHeight
    for j = 1:totWidth
        if sum(IdiffCell{i,j}(:)) ~= 0
            Imotion((i-1)*L+1:i*L,(j-1)*L+1:j*L) = 1;
        end
    end
end
         
figure
colormap gray;
subplot(2,2,1)
imagesc(Iold)
title('Old frame')
axis off;

subplot(2,2,2)
imagesc(Inew)
title('New frame')
axis off;

subplot(2,2,3)
imagesc(Imotion)
title('Motion')
axis off;

subplot(2,2,4)
imagesc(Idiff)
title('Diff')
axis off;

%% ------------------------------------------------------------------------
% -------------------- B TASK 2.4 -----------------------------------------
% -------------------------------------------------------------------------
clc, close all
MV = zeros(totHeight, totWidth, 2);

% Loop through Inew and search in old to find best match
for newPosX=1:totHeight
    for newPosY = 1:totWidth
        if sum(sum(Imotion((newPosX-1)*L+1:newPosX*L,(newPosY-1)*L+1:newPosY*L))) ~= 0
            currBlock=newCell{newPosX,newPosY};
            bestVal=inf;
            for oldPosX = 1:height-16
                for oldPosY = 1:width-16
                    idxX=oldPosX:oldPosX+15;
                    idxY=oldPosY:oldPosY+15;
                    currMAE=sum(sum(abs(currBlock-Iold(idxX,idxY))));
                    if currMAE < bestVal
                        bestVal = currMAE;
                        MV(newPosX,newPosY,1)=oldPosX-(newPosX-1)*L;
                        MV(newPosX,newPosY,2)=oldPosY-(newPosY-1)*L;
                    end
                end
            end
        end
    end
end


Icomp = zeros(height, width);

% I4
I4 = zeros(240,320);
for i=1:totHeight
    for j=1:totWidth
        indxX=((i-1)*L+1:i*L);
        indxY=(j-1)*L+1:j*L;
        if Imotion((i-1)*L+1:i*L,(j-1)*L+1:j*L) == 1
            I4(indxX,indxY)=Iold(indxX+MV(i,j,1),indxY+MV(i,j,2));
        end
    end
end

% I5
I5 = I4;
for i=1:totHeight
    for j=1:totWidth
        indxX=((i-1)*L+1:i*L);
        indxY=(j-1)*L+1:j*L;
        if Imotion((i-1)*L+1:i*L,(j-1)*L+1:j*L) == 0
            I5(indxX,indxY)=Iold(indxX,indxY);
        end
    end
end


error30=30*abs(Inew-I5);

figure
colormap gray
subplot(2,2,1)
% Plot old original fram
imagesc(Inew)
% Plot I4
subplot(2,2,2)
imagesc(I4)
% Plot I5, compressed intra + inter
subplot(2,2,3)
imagesc(I5)
% Plot Error
subplot(2,2,4)
imagesc(error30)

PSNR = psnr(I5, Inew)        
SSIM = ssim(I5, Inew)        
