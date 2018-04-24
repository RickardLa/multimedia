% ------------------------------------------------------------------------
% -------------------- TASK 1.1 -------------------------------------------
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
I2 = rgb2gray(mov(20).frames); 

% imwrite(I1, 'frame10', 'bmp');
% imwrite(I2, 'frame20', 'bmp');

% imagesc(I1), colormap gray
% figure
% imagesc(I2), colormap gray

%% Step 1.2
clc, close

originalImg = mat2gray(I1); 
DCTCoeff = dct2(originalImg); 


% imshow(log(abs(DCTCoeff)),[]), colormap(gca,jet), colorbar


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
imagesc(log(abs(DCTCoeff))), colormap(gca,jet)
title('DCT Domain')
axis off;

saveas(gca, 'Figures/Compressed.eps','epsc');

%% ------------------------------------------------------------------------
% -------------------- TASK 1.2 -------------------------------------------
% ------------------------------------------------------------------------
clc,clear,close;

% Step 2.1
originalImg = mat2gray(imread('Figures/frame10','bmp'));
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
        block = dct2(allBlocks{i,j});
        DCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L) = block; 
    end
end


imagesc(log(abs(DCTBlocks))), colormap(gca,jet)
title('DCT Domain')
axis off;

% 2.3 COMPRESS

th = floor(0.9 * totWidth * totHeight);

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

% Now compute IDCT of all blocks and store them in DCTBlocks
DCTBlocks = zeros(height, width);
for i=1:totHeight       
    for j=1:totWidth   
        block = idct2(allCmprsdBlocks{i,j});
        IDCTBlocks((i-1)*L+1:i*L,(j-1)*L+1:j*L) = block; 
    end
end

errorImg = abs(originalImg - IDCTBlocks);
errorImg30 = 30*abs(originalImg - IDCTBlocks);

PSNR = psnr(IDCTBlocks, originalImg)         % PSNR in dB
SSIM = ssim(IDCTBlocks, originalImg)         % SSIM = 1 means the images are identical

%%
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
imagesc(log(abs(DCTBlocks))), colormap(gca,jet)
title('DCT Domain')
axis off;










