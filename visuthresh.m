%% designed by M.Kiran Kumar
%% R & D , vizag,Andhra Pradesh,INDIA
%% mail: matlabkirankumar@gmail.com


clear all;
close all;
clc;

%Implementing Visu Shrink-
%Denoising using universal threshold with both hard and soft thresholding

%Note: Figure window 1 displays the original image, fig 2 the noisy img
%fig 3 denoised img by hard thresholding, fig 4 denoised by soft thresholding

msgbox('DESIGNED BY M.KIRAN KUMAR: matlabkirankumar@gmail.com')

%Reading the image 
pic=imread('barbara','png');
figure, imagesc(pic);colormap(gray);

%Define the Noise Variance and adding Gaussian noise
%While using 'imnoise' the pixel values(0 to 255) are converted to double in the range 0 to 1
%So variance also has to be suitably converted
sig=100;
V=(sig/256)^2;

npic=imnoise(pic,'gaussian',0,V);
figure, imagesc(npic);colormap(gray);

%Define the type of wavelet(filterbank) used and the number of scales in the wavelet decomp
filtertype='db4';
levels=5;

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


%Soft thresholding
softC=[C(1:S(1,1)^2), sthresh(C(S(1,1)^2+1:length(C)),UT)];

%Reconstructing the image from the soft-thresholded wavelet coefficients
newpics=waverec2(softC,S,filtertype);

%Displaying the soft-denoised image
figure, imagesc(newpics);colormap(gray);
