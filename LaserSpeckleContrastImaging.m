% Author Chang L
% Date: 01-17-19
% Purpose: Laser Speckle Data processing pipeline 
% Data: LSCI measurement of mouse cerebral blood flow during Co2 induced hypercapnia 
% and stroke (carotid artery ligation model); camera exposure time 32 us

%% Get Laser Speckle Contrast Imaging 
clc
clear all
fileName = 'Chang_bl.rls';
startT=0;
T= 0.005; 
Fps =  22800;
totalframe = 98048;  
sizeT=114;  
ROI = [];
type = 'frame';
kernalSize = 7;
i=1;
% Read LSCI datasets
[data2, sampling,timeStamps]=readRLS(fileName,startT,sizeT,ROI,type);
for startT =0: sizeT: totalframe-sizeT  
     [data,sampling,timeStamps] = readRLS(fileName,startT,sizeT,ROI,type);
     speckle(:,:,i) = squeeze(sum(data,3));
      i=i+1
end
% Get spatial contrast images
 for i=1:size(speckle,3)
  LSCIimg(:,:,i) =  SLSCI (speckle(:,:,i),kernalSize);
  i=i+1
 end
LSCIavg= squeeze(mean(LSCLimg,3)); 
flowmap = 1 ./ (LSCIavg.*LSCIavg) ;
imagesc(flowmap)
colorbar
colormap jet
axis image
caxis([prctile(flowmap(:),5) prctile(flowmap(:),95)]);
title('baseline flowmap');
%Data preprocessing crop the LSCI image if it has 960 X 1024 pixels.
rect = [1,371,1024,149];
LSCIblroi=imcrop(LSCIbl,rect); 
%% Mask image to identify different vessel types.
[finalMask, meanMap] = Mask (flowmap_cropped);
% Get CBF with the conventional 1/K^2 model.
lv=meanMap(find(meanMap.* (finalMask==5)));
CBFlv =mean(lv);      SDlv=std(lv);
sv=meanMap(find(meanMap.*(finalMask==3)));
CBFsv =mean(sv);      SDsv=std(sv);
par=meanMap(find(meanMap.*(finalMask==1)));
CBFpar =mean(par);    SDpar=std(par);
%%
--------------------------------Self-defined functions -----------
% SLSCI: spatial contrast image
function [sLSCI] = SLSCI(Imagedata,kernalSize)

Kernal=single((ones(kernalSize, kernalSize)));
frameStd = stdfilt(Imagedata,Kernal); 
frameMean= conv2(single(Imagedata),(Kernal),'same')  /  sum(Kernal(:));
sLSCI=frameStd./frameMean;

end
% Mask to classify different vessel types
function [finalMask, meanMap] = Mask (map)
meanMap=map(:,200:end);
typesMask=zeros(size(meanMap));
typesMask(meanMap>prctile(meanMap(:),5))=1;  
typesMask(meanMap>prctile(meanMap(:),54))=2;  % 1 is the parenchyma:
typesMask(meanMap>prctile(meanMap(:),60))=3;    
typesMask(meanMap>prctile(meanMap(:),84))=4;
typesMask(meanMap>prctile(meanMap(:),90))=5;
typesMask(meanMap>prctile(meanMap(:),99.5))=6;    % 94.98%   65.22%
figure;imagesc(typesMask); 
finalMask=zeros(size(typesMask)); 

disp('Manual mask for the parenchyma')
answ=1;
while answ==1
[x,y]=ginput(4);
BW = poly2mask(x,y,size(finalMask,1),size(finalMask,2));
finalMask(BW==1 & typesMask==1)=1;
answ=input('One more region for parenchyma?');
end

disp('Manual mask for the small vessels')
answ=1;
while answ==1
[x,y]=ginput(4);
BW = poly2mask(x,y,size(finalMask,1),size(finalMask,2));
finalMask(BW==1 & typesMask==3)=3;
answ=input('One more region for the small vessels?');
end

disp('Manual mask for the large vessels')
answ=1;
while answ==1
[x,y]=ginput(4);
BW = poly2mask(x,y,size(finalMask,1),size(finalMask,2));
finalMask(BW==1 & typesMask==5)=5;
answ=input('One more region for the large vessels?');
end

end

% readRLS from BUNPC/laserSpeckleImaging/Raw/readRLS.m --written by ddpostnov
function [data,sampling,timeStamps]=readRLS(fileName,startT,sizeT,ROI,type)
fileReadId = fopen(fileName, 'r');  % file identifier
fseek(fileReadId,0*1024,-1 );  % -1: read from the beginning of file
%NOTE IN SOME RLS FILES Y and X were misplaced!!! Check the image and fix
%it!
% ROI = [xmin  width
%        ymin  height]
sizeX=fread(fileReadId,1,'*uint64');
sizeY=fread(fileReadId,1,'*uint64');
 % note that I exchange the value of sizeX and sizeY 
if sizeT==0
    sizeT=fread(fileReadId,1,'*uint64'); % if sizeT=0, then read all data. 98048 frames
else
    fread(fileReadId,1,'*uint64'); 
end
sampling=fread(fileReadId,1,'*uint64');
timeStamps=zeros(sizeT,1,'int64');
firstByte=30*1024+sizeX*sizeY*startT+8*startT; % start T ~ timestamps
% T should start at 0, so that the first frame is readed. if startT=1,
% it's skipped! 
fseek(fileReadId,firstByte,-1 );  
if strcmp(type,'frame')
    data=zeros(sizeX,sizeY,sizeT,'uint8');
    for t=1:1:sizeT
        timeStamps(t)=fread(fileReadId,1,'*uint64');
        data(:,:,t)=fread(fileReadId,[sizeX,sizeY],'*uint8');
    end
elseif strcmp(type,'rect')
    data=zeros(length(ROI(1,1):1:ROI(1,2)),length(ROI(2,1):1:ROI(2,2)),sizeT,'uint8');
    for t=1:1:sizeT
        timeStamps=fread(fileReadId,1,'*uint64');
        frame=fread(fileReadId,[sizeX,sizeY],'*uint8');
        data(:,:,t)=frame(ROI(1,1):1:ROI(1,2),ROI(2,1):1:ROI(2,2));
    end
elseif strcmp(type,'1d')
    data=zeros(length(ROI),sizeT,'uint8');
%        data = zeros(sizeT,1);
    for t=1:1:sizeT
        timeStamps=fread(fileReadId,1,'*uint64');
        frame=fread(fileReadId,[sizeX,sizeY],'*uint8');
       data(:,t)=frame(ROI);
% %         frame = imcrop(frame,rect = [1,371,1024,150]);  % only for baseline
% %        data(t,1)= frame(row,col);
    end
elseif strcmp(type,'ROI')
    row = ROI(:,1);
    col = ROI(:,2);
    data =zeros(length(ROI(:,1)),sizeT,'uint8');
    for t=1:1:sizeT
        t
        timeStamps =fread(fileReadId,1,'*uint64');
        frame =fread(fileReadId,[sizeX,sizeY],'*uint8');
        for i=1:length(row)
        data(i,t)=frame(row(i),col(i));
        end
    end   
end
fclose(fileReadId);

end














