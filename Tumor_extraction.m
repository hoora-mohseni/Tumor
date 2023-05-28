function Tumor_extraction


% In this project we propose a simple-to-implement and fast strategy for 
% detection and extraction of brain tumors in Magnetic Resonance Imaging
% (MRI) of the brain. The proposed method consists of two steps. The first 
% step is preprocessing of given MRI image using high-pass, low-pass and 
% median filters and in the second step we use methods such as 
% thresholding, watershed algorithm, and morphological operation. In the 
% end our output? will be a tumor region that we can calculate its area and 
% central point using morphological operations
%
% **********************************************************************
%
% Running this code will want you to :
%
% 1) select a Protocol number which is the extraction technique discussed 
% in the report.
% 2) input image, which you want to extract the tumor area from it.
%
% **********************************************************************


%% Input

% Selecting the Protocol Number.

Font_size = get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 14);
Protocol_num = menu('Choose a Protocol Number', '1', '2', '3');
set(0, 'DefaultUIControlFontSize', Font_size);

% Load the input file.

[I,path] = uigetfile('*.jpg','Select an input image');
str      = strcat(path,I);
img      = imread(str);

%% Section 1

% Calculating run-time of the code.

tic 

% Convert the image into grayscale.

img = rgb2gray(img);

% Convert the image into double format.

img = im2double(img);

% Resize the image so that all the images have the same size.

img = imresize(img,[630 500]);

% Normalize all the images in [0,1] intensity values.

img = img - min(img(:));
img = img ./ max(img(:));

[m,n] = size(img);

% Applying a high-pass Laplacian filter.

h              = fspecial('laplacian');
Laplacian_mask = imfilter(img,h);
img2           = img + (-1).*Laplacian_mask;

% Applying a Median fitler of size 3*3

img3 = medfilt2(img2,[3 3]);

% Normalizing the resultant image so to binarize it in the next step.

img4 = img3;
img4 = img4 - min(img4(:));
img4 = img4 ./ max(img4(:));

%% Section 2

% Binarizing the resultant image.

thresh = mean(img4(:));
img5 = im2bw(img4, 1.75*thresh);

switch Protocol_num
    
        case 1

        % Defining Structuring Element
        se   = strel('disk',15);

        % Erosion
        img5_erode    = imerode(img5,se);
        
        % Dilation
        img5_dilate   = imdilate(img5_erode,se);
        
    case 2
        
        % Defining Structuring Element
        se   = strel('disk',15);
        
        % Opening by reconstruction
        img5_erode    = imerode(img5,se);
        img5_open_rec = imreconstruct(img5_erode,img5);

        % Closing by reconstruction
        img5_dilate    = imdilate(img5_open_rec,se);
        img5_close_rec = imreconstruct(imcomplement(img5_dilate),imcomplement(img5_open_rec));
        img5_close_rec = imcomplement(img5_close_rec);

    case 3
        
        
        % Extracting Connected Components
        CC = bwconncomp(img5);
        numOfPixels = cellfun(@numel,CC.PixelIdxList);

        % Calculating Area and Perimeter of each connected component
        stats_Area      = regionprops(CC,'Area');
        stats_Perimeter = regionprops(CC,'Perimeter');

        % Calculating Circularity of each connnected component
        Circularity = ([stats_Area.Area].*4.*pi) ./ ([stats_Perimeter.Perimeter].*[stats_Perimeter.Perimeter]+eps);

        % Excluding the components with numOfPixels under 5000 and 
        % Circularity over 0.1.
        for i = 1 : length(numOfPixels)
            if numOfPixels(i) >= 5000 && Circularity(i) < 0.1
                numOfPixels_new(i) = 0;
            else
                numOfPixels_new(i) = numOfPixels(i);
            end
        end

        % Extracting the largest connected component which meets two
        % conditions defined above.
        [~,indexOfMax] = max(numOfPixels_new);
        biggest = zeros(size(img5));
        biggest(CC.PixelIdxList{indexOfMax}) = 1;

end
%% Plots

% Plotting the resultant images.

figure;
subplot(1,2,1);
imshow(img);
title('Original Image');
subplot(1,2,2);
imshow(img2);
title('Image + Laplacian mask');

figure;
subplot(1,2,1);
imshow(img3);
title('Image + Laplacian mask + Median filter');
subplot(1,2,2);
imshow(img5);
title('Binarized Image');

switch Protocol_num
    
    case 1
        
        figure;
        subplot(1,2,1);
        imshow(img5_erode);
        title('Erosion');
        subplot(1,2,2);
        imshow(img5_dilate);
        title('Eroded Image followed by Dilation');

    case 2

        figure;
        subplot(1,2,1);
        imshow(img5_open_rec);
        title('Opening by reconstruction');
        subplot(1,2,2);
        imshow(img5_close_rec);
        title('Closing by reconstruction');
        
    case 3
        
        figure;
        subplot(1,2,1);
        imshow(img5);
        title('Binarized Image');
        subplot(1,2,2);
        imshow(biggest);
        title('Tumor Extracted');
        
end

%% Boundary Extraction + Plot

% Extract the boundary of the tumor region.

switch Protocol_num
    case 1
        
        A    = bwboundaries(img5_dilate);
        B    = A{1};
        img6 = imfuse(img,img5_dilate);
        img6 = rgb2gray(img6);
        img6 = im2double(img6);
        
    case 2
        
        A    = bwboundaries(img5_close_rec);
        B    = A{1};
        img6 = imfuse(img,img5_close_rec);
        img6 = rgb2gray(img6);
        img6 = im2double(img6);
        
    case 3
        
        A = bwboundaries(biggest);
        B    = A{1};
        img6 = imfuse(img,biggest);
        img6 = rgb2gray(img6);
        img6 = im2double(img6);
        
end



% Plotting the boundary of the tumor region superimposed on original image.

figure;
subplot(1,2,1);
imshow(img);
title('Original image');
subplot(1,2,2);
imshow(img6);
title('Bright tumor area superimposed on the original image');
hold on
plot(B(:,2),B(:,1),'c:','Linewidth',2);

%% Calculation of Area and Centroid of each tumor + Plot

% Calculating Area and Centroid of the tumor regions using regionprops.

switch Protocol_num
    case 1
        
        s = regionprops(img5_dilate,'area');
        c = regionprops(img5_dilate,'centroid');
        
    case 2
        
        s = regionprops(img5_close_rec,'area');
        c = regionprops(img5_close_rec,'centroid');
        
    case 3
        
        s = regionprops(biggest,'area');
        c = regionprops(biggest,'centroid');
        
end

% Superimposing Centroid of the tumor region and it's Area on original
% image.

hold on
plot(c.Centroid(1),c.Centroid(2),'g+','Markersize',25);
message = sprintf('The area of tumor is equal to %0.2f %% of the image',100* s.Area/(m*n));
text(c.Centroid(1)-350,c.Centroid(2)+100,message,'color','red','Fontsize',14,'Fontweight','Bold');

toc

end