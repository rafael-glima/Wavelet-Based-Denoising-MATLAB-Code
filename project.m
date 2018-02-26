%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    EE530/IMAGE PROCESSING - FINAL PROJECT             %
%                        By: Rafael Lima - ID: 20155403                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  This is a matlab implementation of the wavelet-based denoising algorithm 
%  proposed in the paper entitled "Wavelet Image Threshold Denoising Based
%  on Edge Detection", by Wei Liu and Zhengming Ma (IMACS 2006). It was tes-
%  ted with a 256x256 "Lena" image. An auxiliary code with a Visushrink
%  threshold algorithm implementation, written by Kiran Kumar and posted in
%  the Mathworks website (www.mathworks.com) was also used in the end of 
%  this file for comparisons.

clear all
close all
clc

% Loading image and getting the Y component

Img = imread('lena.jpg');
Img = rgb2gray(Img);

% Defining the strength of the noise

sig = 10;
noise_var = (sig/256)^2;

%%%%%%%%% Adding Noise %%%%%%%

Noise_Img = imnoise(Img,'gaussian',0,noise_var);


% Plotting original image and the the noisy version

figure
imagesc(Img);colormap(gray);
title('Original Image')

pause(0.01)

figure
imagesc(Noise_Img);colormap(gray);
title('Image with noise')

%%%%%%%%%% Averaging Filter %%%%%%%%%%%%

%Noise_Img = imfilter(Noise_Img,ones(2)/4,'symmetric','same');

%%%%%%%%%% Wavelet Coefficients of each row and collumn %%%%%%%%%%%

Noise_Img = double(Noise_Img);
s = size(Noise_Img);

H = 1.4142*[0.025 0.475 0.475 0.025];
G = 1.4142*[0.025 -0.475 0.475 -0.025];

beta = 0.2; % Threshold parameter


% 2-D Wavelet decomposition of the image (since the number of levels used in
% the denoising method was not specified, the maximum number of levels was
% used). Also, the wavelet used was the Daubechie 2 ('db2').

levels = log2(size(Noise_Img,1)); 

[C,S] = wavedec2(Noise_Img,levels,'db2')%H,G);

lip_limit = 0; % Limit of lipschitz exponents for thresholding

% Extraction and processing of the components in the vector corresponding
% to the horizontal decompositions. This "for" loop was implemented for
% recursively getting the indexes of beginning and ending of the total 8
% levels of horizontal decomposition components and, at the same time,
% calculating the lipschitz exponents of each of its elements, classifying
% them as edges or noise, and applying distinct thresholds for each of
% these two categories (soft-threshold for edge-related components and hard
% -threshold for remaining ones).

start_ind = 0;

for i = 1:levels
    i
    
    % Calculating the threshold according to the level
    
    thresh = sig*sqrt(4*log(prod(S(i+1,:))));
    
   if i == 1
        start_ind = start_ind + prod(S(i,:));
        end_ind = start_ind + prod(S(i+1,:));  
        H_Coeffs = C(start_ind+1:end_ind);
        C(start_ind+1:end_ind) = H_Coeffs;
   elseif i == 8
        start_ind = start_ind + 3*prod(S(i,:));
        end_ind = start_ind + prod(S(i+1,:));
        H_Coeffs = reshape(H_Coeffs,sqrt(length(H_Coeffs)),sqrt(length(H_Coeffs)));
        H_Coeffs = H_Coeffs(2:end-1,2:end-1);
        
        % Storing the coefficients of the adjacent level, for calculating
        % the Lipschitz exponents according to the description in equation (9) of 
        % the paper entitled "Measurement of Lipschitz Exponent (LE) using Wavelet
        % Transform Modulus Maxima (WTMM)", by Venkatakrishnan et al.
        % (IJSER - June/2012).
        
        next_lvl_H_Coeffs = kron(H_Coeffs,ones(2));
        next_lvl_H_Coeffs = padarray(next_lvl_H_Coeffs,[1 1],'symmetric','post');
        next_lvl_H_Coeffs = reshape(next_lvl_H_Coeffs,1,size(next_lvl_H_Coeffs,1)*size(next_lvl_H_Coeffs,2));
        H_Coeffs = C(start_ind+1:end_ind);
        lipschitz_exp = log2(abs(next_lvl_H_Coeffs./H_Coeffs));
   else
        start_ind = start_ind + 3*prod(S(i,:));
        end_ind = start_ind + prod(S(i+1,:));
        H_Coeffs = reshape(H_Coeffs,sqrt(length(H_Coeffs)),sqrt(length(H_Coeffs)));
        H_Coeffs = H_Coeffs(2:end-1,2:end-1);
        next_lvl_H_Coeffs = kron(H_Coeffs,ones(2));
        next_lvl_H_Coeffs = padarray(next_lvl_H_Coeffs,[1 1],'symmetric');
        next_lvl_H_Coeffs = reshape(next_lvl_H_Coeffs,1,size(next_lvl_H_Coeffs,1)*size(next_lvl_H_Coeffs,2));
        H_Coeffs = C(start_ind+1:end_ind);
        lipschitz_exp = log2(abs(next_lvl_H_Coeffs./H_Coeffs));
   end
   
   % Classification of the components according to the value of their
   % Lipschitz Exponents ( lipschitz_exp > 0 : edge ; lipschitz_exp <= 0:
   % noisy component).
   
   if i > 1
         for j = 1:length(H_Coeffs)
            if (lipschitz_exp(j) > lip_limit) && (abs(H_Coeffs(j)) < beta*thresh)
                H_Coeffs(j) = 0;
            elseif (lipschitz_exp(j) <= lip_limit) && (abs(H_Coeffs(j)) < thresh)
               H_Coeffs(j) = 0;
           end
         end
   else
       for j= 1:length(H_Coeffs)
            if (abs(H_Coeffs(j)) < thresh)
                H_Coeffs(j) = 0;
            end
       end
   end
   
   C(start_ind+1:end_ind) = H_Coeffs;
   
end

% Here, the same procedure of separation, classification and categorical 
% thresholding is now applied to the components related to vertical
% decomposition.

start_ind = 0;

for i = 1:levels
    i
    thresh = sig*sqrt(4*log(prod(S(i+1,:))));
   if i == 1
        start_ind = start_ind + prod(S(i,:));
        end_ind = start_ind + prod(S(i+1,:));  
        V_Coeffs = C(start_ind+1:end_ind);
        thresh = 0;
        C(start_ind+1:end_ind) = V_Coeffs;
   elseif i == 8
        start_ind = start_ind + 3*prod(S(i,:));
        end_ind = start_ind + prod(S(i+1,:));
        V_Coeffs = reshape(V_Coeffs,sqrt(length(V_Coeffs)),sqrt(length(V_Coeffs)));
        V_Coeffs = V_Coeffs(2:end-1,2:end-1);
        next_lvl_V_Coeffs = kron(V_Coeffs,ones(2));
        next_lvl_V_Coeffs = padarray(next_lvl_V_Coeffs,[1 1],'symmetric','post');
        next_lvl_V_Coeffs = reshape(next_lvl_V_Coeffs,1,size(next_lvl_V_Coeffs,1)*size(next_lvl_V_Coeffs,2));
        V_Coeffs = C(start_ind+1:end_ind);
        lipschitz_exp = log2(abs(next_lvl_V_Coeffs./V_Coeffs));
   else
        start_ind = start_ind + 3*prod(S(i,:));
        end_ind = start_ind + prod(S(i+1,:));
        V_Coeffs = reshape(V_Coeffs,sqrt(length(V_Coeffs)),sqrt(length(V_Coeffs)));
        V_Coeffs = V_Coeffs(2:end-1,2:end-1);
        next_lvl_V_Coeffs = kron(V_Coeffs,ones(2));
        next_lvl_V_Coeffs = padarray(next_lvl_V_Coeffs,[1 1],'symmetric');
        next_lvl_V_Coeffs = reshape(next_lvl_V_Coeffs,1,size(next_lvl_V_Coeffs,1)*size(next_lvl_V_Coeffs,2));
        V_Coeffs = C(start_ind+1:end_ind);
        lipschitz_exp = log2(abs(next_lvl_V_Coeffs./V_Coeffs));
   end
   
   if i > 1
         for j = 1:length(V_Coeffs)
            if (lipschitz_exp(j) > lip_limit) && (abs(V_Coeffs(j)) < beta*thresh)
                V_Coeffs(j) = 0;
            elseif (lipschitz_exp(j) <= lip_limit) && (abs(V_Coeffs(j)) < thresh)
               V_Coeffs(j) = 0;
           end
         end
   else
       for j= 1:length(V_Coeffs)
            if (abs(V_Coeffs(j)) < thresh)
                V_Coeffs(j) = 0;
            end
       end
   end
   
   C(start_ind+1:end_ind) = V_Coeffs;
   
end

% % Here, the same procedure of separation, classification and categorical 
% thresholding is now applied to the components related to diagonal
% decomposition.

start_ind = 0;

for i = 1:levels
    i
    thresh = sig*sqrt(4*log(prod(S(i+1,:))));
   if i == 1
        start_ind = start_ind + prod(S(i,:));
        end_ind = start_ind + prod(S(i+1,:));  
        D_Coeffs = C(start_ind+1:end_ind);
        thresh = 0;
        C(start_ind+1:end_ind) = D_Coeffs;
   elseif i == 8
        start_ind = start_ind + 3*prod(S(i,:));
        end_ind = start_ind + prod(S(i+1,:));
        D_Coeffs = reshape(D_Coeffs,sqrt(length(D_Coeffs)),sqrt(length(D_Coeffs)));
        D_Coeffs = D_Coeffs(2:end-1,2:end-1);
        next_lvl_D_Coeffs = kron(D_Coeffs,ones(2));
        next_lvl_D_Coeffs = padarray(next_lvl_D_Coeffs,[1 1],'symmetric','post');
        next_lvl_D_Coeffs = reshape(next_lvl_D_Coeffs,1,size(next_lvl_D_Coeffs,1)*size(next_lvl_D_Coeffs,2));
        D_Coeffs = C(start_ind+1:end_ind);
        lipschitz_exp = log2(abs(next_lvl_D_Coeffs./D_Coeffs));
   else
        start_ind = start_ind + 3*prod(S(i,:));
        end_ind = start_ind + prod(S(i+1,:));
        D_Coeffs = reshape(D_Coeffs,sqrt(length(D_Coeffs)),sqrt(length(D_Coeffs)));
        D_Coeffs = D_Coeffs(2:end-1,2:end-1);
        next_lvl_D_Coeffs = kron(D_Coeffs,ones(2));
        next_lvl_D_Coeffs = padarray(next_lvl_D_Coeffs,[1 1],'symmetric');
        next_lvl_D_Coeffs = reshape(next_lvl_D_Coeffs,1,size(next_lvl_D_Coeffs,1)*size(next_lvl_D_Coeffs,2));
        D_Coeffs = C(start_ind+1:end_ind);
        lipschitz_exp = log2(abs(next_lvl_D_Coeffs./D_Coeffs));
   end
   
   if i > 1
         for j = 1:length(D_Coeffs)
            if (lipschitz_exp(j) > lip_limit) && (abs(D_Coeffs(j)) < beta*thresh)
                D_Coeffs(j) = 0;
            elseif (lipschitz_exp(j) <= lip_limit) && (abs(D_Coeffs(j)) < thresh)
               D_Coeffs(j) = 0;
           end
         end
   else
       for j= 1:length(D_Coeffs)
            if (abs(D_Coeffs(j)) < thresh)
                D_Coeffs(j) = 0;
            end
       end
   end
   
   C(start_ind+1:end_ind) = D_Coeffs;
   
end

% Reconstructing Image with the thresholded components

Rec_Img = waverec2(C,S,'db2')%H,G);
Rec_Img = uint8(Rec_Img);

figure
imagesc(Rec_Img);colormap(gray);
title('Recovered Denoised Image')
pause(0.01)

% Calculating and printing the SNR of both the noisy and the denoised
% picture

[peaksnr, snr] = psnr(uint8(Noise_Img),Img);

fprintf('\n The peak-SNR value of Noisy Image is %0.4f', peaksnr);
fprintf('\n The SNR value of Noisy Image is %0.4f', snr);

[peaksnr, snr] = psnr(Rec_Img,Img);

fprintf('\n The peak-SNR value of Recovered Image is %0.4f', peaksnr);
fprintf('\n The SNR value of Recovered Image is %0.4f \n', snr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Visushrink threshold: This portion of the code is NOT MINE. It was   %
%    downloaded from the Mathworks website for comparison through the     %
%    use of the same exact values for beta (threshold parameter) and      %
%    sigma (noise strength).                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pic = Img;

%figure, imagesc(pic);colormap(gray);

%Define the Noise Variance and adding Gaussian noise
%While using 'imnoise' the pixel values(0 to 255) are converted to double in the range 0 to 1
%So variance also has to be suitably converted

V=(sig/256)^2;

npic = uint8(Noise_Img);

%npic=imnoise(pic,'gaussian',0,V);
%figure, imagesc(npic);colormap(gray);

%Define the type of wavelet(filterbank) used and the number of scales in the wavelet decomp
filtertype='db2';
%levels=5;

%Doing the wavelet decomposition
[C,S]=wavedec2(npic,levels,filtertype);

%Define the threshold(universal threshold)
M=size(pic,1)^2;
UT=sig*sqrt(2*log(M));

%Hard thresholding
%Doing the hard thresholding-threshold only detail coefficients!!
hardC=[C(1:S(1,1)^2), hthresh(C(S(1,1)^2+1:length(C)),UT)];

%Reconstructing the image from the hard-thresholded wavelet coefficients
newpich=waverec2(hardC,S,filtertype);

%Displaying the hard-denoised image
figure, imagesc(newpich);colormap(gray);
title('Hard-thresholded Visushrink Recovered Image')
pause(0.01)


%Soft thresholding
softC=[C(1:S(1,1)^2), sthresh(C(S(1,1)^2+1:length(C)),UT)];

%Reconstructing the image from the soft-thresholded wavelet coefficients
newpics=waverec2(softC,S,filtertype);

%Displaying the soft-denoised image
figure, imagesc(newpics);colormap(gray);
title('Soft-thresholded Visushrink Recovered Image')
pause(0.01)

% Calculating and printing the SNR values of the two denoised pictures
% (hard-threshold and soft-threshold)

[peaksnr, snr] = psnr(uint8(newpich),pic);

fprintf('\n The peak-SNR value of Hard-Thresholded Recovered Image is %0.4f', peaksnr);
fprintf('\n The SNR value of Hard-Thresholded Recovered Image is %0.4f', snr);

[peaksnr, snr] = psnr(uint8(newpics),pic);

fprintf('\n The peak-SNR value of Soft-Thresholded Recovered Image is %0.4f', peaksnr);
fprintf('\n The SNR value of Soft-Thresholded Recovered Image is %0.4f \n', snr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
