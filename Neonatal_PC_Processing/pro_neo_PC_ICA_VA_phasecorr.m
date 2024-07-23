% script for PC MRI processing
% calculate the fux (ml/min) of internal cortid a. and vetebra a.
% hand drawing blood vessel mask based on the magnitude(for vessels) image

% Shin-Lei Peng
% 03/25/2015
% modified by Dengrong Jiang 07/12/2018: extract FOV and VENC from
% dicominfo of phase images instead of hard-coding

close all;
clear;
% Quickmode: in this mode, the user only needs to choose the folder of 
% "MAG" image, the folder of the "P" image will be computed automatically.
% However, to use this mode, the image files need to be organized in a
% certain way, as if they are exported from the scanner to the root
% directory of a thumb drive or hard driven
bQuickMode = 0;

% magnitude image
if bQuickMode
    % To save time, but requires a specific organization of images
    pathNamePar = uigetdir('Choose the "MAG" folder:');
    if pathNamePar(end) ~= filesep
        pathNamePar = [pathNamePar filesep];
    end
    P = dir(fullfile(pathNamePar,'*.ima'));
    if isempty(P)
        P = dir(fullfile(pathNamePar,'*.dcm'));
    end
    slicenum=size(P,1);
    for i=1:slicenum
        magname=P(i,:).name
        h=dicominfo(fullfile(pathNamePar, magname));
        tmp=dicomread(h);
        pcimg(:,:,h.InstanceNumber)=tmp;
    end
else
    [fileNamePar,pathNamePar] = uigetfile({'*.ima';'*.dcm'},'Choose the dicom image in the "MAG" folder:','MultiSelect','off');
    [~,~,fname_ext] = fileparts([pathNamePar fileNamePar]);
    h=dicominfo([pathNamePar, fileNamePar]);
    pcimg=dicomread(h);
end

newpcimg=RealignImageArray(pcimg,'-y+x+z');  % rl, ap, and fh directions
matrix=[size(newpcimg), 1];
resolution=[h.PixelSpacing' h.SliceThickness];
write_ANALYZE(newpcimg,[pathNamePar 'magnitude.img'],matrix, resolution);
vesimage=double(pcimg); % display in DICOM orientation

%phase image
if bQuickMode
    % To save time, but requires a specific organization of images
    mag_num = str2num(pathNamePar(end-4:end-1));
    if mag_num == 98
        pha_num = 100; % number 99 is always PHOENIXZIPREPORT_0099
    else
        pha_num = mag_num + 1;
    end
    if pha_num < 10
        pha_num_str = ['000', num2str(pha_num)];
    elseif pha_num < 100
        pha_num_str = ['00', num2str(pha_num)];
    elseif pha_num < 1000
        pha_num_str = ['0', num2str(pha_num)];
    elseif pha_num < 10000
        pha_num_str = num2str(pha_num);
    else
        fprintf('Image series number is out of range!\n');
        return;
    end
    pathNamePar = [pathNamePar(1:end-9), 'P_', pha_num_str, filesep];
    P = dir(fullfile(pathNamePar,'*.ima'));
    if isempty(P)
        P = dir(fullfile(pathNamePar,'*.dcm'));
    end
    [fname_path,fname_body,fname_ext] = fileparts([pathNamePar, filesep, P.name]);
    filename = fullfile(fname_path,fname_body);
    slicenum=size(P,1);
    for i=1:slicenum
        phasename=P(i,:).name
        h=dicominfo(fullfile(pathNamePar, phasename));
        tmp=dicomread(h);
        pcimg(:,:,h.InstanceNumber)=tmp;
    end
else
    [fileNamePar,pathNamePar] = uigetfile({'*.ima';'*.dcm'},'Choose the dicom image in the "P" folder:','MultiSelect','off');
    [fname_path,fname_body,fname_ext] = fileparts([pathNamePar fileNamePar]);
    filename = fullfile(fname_path,fname_body);
    h=dicominfo([pathNamePar, fileNamePar]);
    pcimg = dicomread(h);
end

% get the VENC and FOV from dicominfo
venc_label = h.Private_0051_1014;
number_end = strfind(venc_label, '_');
Venc = str2double(venc_label(2:number_end(1)-1));
fov_label = h.Private_0051_100c;
times_sign = strfind(fov_label, '*');
FOV = zeros(1,2);
FOV(1) = str2double(fov_label(4:times_sign(1)-1));
FOV(2) = str2double(fov_label(times_sign(1)+1:end));
% verify VENC and FOV
answer = questdlg(sprintf('VENC=%d, FOV=%d*%d', Venc, FOV(1), FOV(2)), ...
	'Are these correct?', ...
	'Yes','No','Yes');
if strcmp(answer, 'No') % if VENC and FOV are wrong, let user input
    Venc = str2double(inputdlg('VENC (cm/s)'));
    fov_answer = inputdlg('FOV (mm), e.g. 100*100');
    fov_label = fov_answer{1};
    times_sign = strfind(fov_label, '*');
    FOV = zeros(1,2);
    FOV(1) = str2double(fov_label(1:times_sign(1)-1));
    FOV(2) = str2double(fov_label(times_sign(1)+1:end));
end

newpcimg=RealignImageArray(pcimg,'-y+x+z');  % rl, ap, and fh directions
matrix=[size(newpcimg),1];
resolution=[h.PixelSpacing' h.SliceThickness];
write_ANALYZE(newpcimg,[pathNamePar 'phaseimg.img'],matrix, resolution);
phaseimage=double(pcimg)*2*Venc/4096-Venc; % display in DICOM orientation
write_ANALYZE(RealignImageArray(phaseimage,'-y+x+z'),[pathNamePar 'velocityimg.img'],matrix, resolution);

[row, col]=size(phaseimage);

index=input('artery name:  ','s');

%%%%%%%%%%%%%%%%%%%%%%%
% ROI drawing
%%%%%%%%%%%%%%%%%%%%%%%
bkg=vesimage(142:163, 100:121);% 250:301,51:100);%
bkg=mean(bkg(:))
show_imgs_sc(vesimage,1,0, 3*bkg);impixelinfo; % for "Sum of squares"
daspect([h.PixelSpacing; 1])
addToolbarExplorationButtons(gcf);
title('Draw Artery ROI');
[mask,xi,yi] = roipoly(); % draw mask
hold on
plot(xi,yi,'r')
hold off
saveas(gcf, [pathNamePar, index, '_ROI.fig']);

%%%%%%%%%%%%%%%%%%%%%%%
% Phase unwrapping
%%%%%%%%%%%%%%%%%%%%%%%
show_imgs_sc(phaseimage,1,-10,10);impixelinfo;
hold on; plot(xi,yi,'r'); hold off;
title('Phase wrapping inside artery? Go to command window');
daspect([h.PixelSpacing; 1])
addToolbarExplorationButtons(gcf);

phasecorr=input('Phase unwrapping? (y/n):  ','s');
if phasecorr == 'y'
    mask_phase_wrap = roipoly();
    cphaseimage=phaseimage;
    for i=1:row
        for j=1:col
            if mask_phase_wrap(i,j)==1
                if phaseimage(i,j)<-0.1
                    cphaseimage(i,j)=phaseimage(i,j)+2*Venc;                
                end
            end
        end
    end
    write_ANALYZE(RealignImageArray(cphaseimage,'-y+x+z'),[fname_path filesep 'PCAP_corr'],matrix,resolution,1,16);
    phaseimage=cphaseimage;    
end


%%%%%%%%%%%%%%%%%%%%%%%
% Compute flux
%%%%%%%%%%%%%%%%%%%%%%%
pixnum_yn=0;% include nagative pixel in phaseimage
flux_yn=0;% include nagative pixel in phaseimage
v_yn=0;% include nagative pixel in phaseimage
mag_yn=0;
vessel_yn=zeros(row, col); % the final vessel area, include nagative pixel in phaseimage
for i=1:row
    for j=1:col
        if mask(i,j)==1 && phaseimage(i,j) > 0
            vessel_yn(i,j)=1;
            pixnum_yn=pixnum_yn+1;
            flux_yn=flux_yn+phaseimage(i,j);
            v_yn=[v_yn;phaseimage(i,j)];
        end
    end
end
show_imgs_sc(phaseimage,1,-10,10);
daspect([h.PixelSpacing; 1])
colormap(gray)
hold on
plot(xi,yi,'r')
axis off
axis square
saveas(gcf, [pathNamePar, index, '_ROI_velocity.fig']);

%------------------------------------calculation
meanv_yn=flux_yn/pixnum_yn; % include negative part
v_yn=v_yn(2:length(v_yn));
maxv_yn=max(v_yn);
flux_yn=flux_yn*60*FOV(1)*FOV(2)/row/col/100;
area_yn=pixnum_yn*FOV(1)*FOV(2)/row/col/100;

%-------------------------------------display
% close all;
fprintf('include negative pixels in phaseimage\n');
flux_yn
area_yn
pixnum_yn
meanv_yn
maxv_yn

output_array = [bkg, flux_yn, area_yn, pixnum_yn, meanv_yn, maxv_yn, Venc, h.SliceThickness]; % more convenient to copy to excel

%% Write ROI results into .csv file
fresult=[fname_path filesep 'flux_' fname_body '_' index '.csv'];
T = table(bkg, flux_yn, area_yn, pixnum_yn, meanv_yn, maxv_yn, Venc, h.SliceThickness,...
    'VariableNames', ["background pixel intensity", "flux (ml/min)", "area (cm^2)", "number of pixels",...
    "mean velocity (cm/s)", "max velocity (cm/s)", "VENC (cm/s)", "slice (mm)"]);
writetable(T, fresult, 'WriteRowNames', false, 'WriteVariableNames', true);

%% Write .txt
    fresult=[fname_path filesep 'flux_' fname_body '_' index '.txt'];
    f=fopen(fresult, 'wt');
    fprintf(f,'%s',filename); 
    fprintf(f,'\n');
    fprintf(f,'no threshold');
    fprintf(f,'\n');
    fprintf(f,'Include negative voxels:');
    fprintf(f,'\n');
      
        fprintf(f,'%s',index);
        fprintf(f,' flux (ml/min): ');
        fprintf(f,'%1.4f',flux_yn);
        fprintf(f,'\n');   
        fprintf(f,'area: ');
        fprintf(f,'%1.4f',area_yn);
        fprintf(f,'\n');   
        fprintf(f,'pixel#: ');
        fprintf(f,'%3.0f',pixnum_yn);
        fprintf(f,'\n');   
        fprintf(f,'mean_v(cm/s): ');
        fprintf(f,'%1.4f',meanv_yn);
        fprintf(f,'\n');   
        fprintf(f,'max_v(cm/s): ');
        fprintf(f,'%1.4f',maxv_yn);
        fprintf(f,'\n');   
    
     fclose(f);
