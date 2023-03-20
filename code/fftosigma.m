% This is the main MATLAB processing file used to get 3D "spectroscopic
% images" (1st and 2nd dims: x,y axes; 3rd dim: k)

% The 2D Fourier variables are the direction cosines d1 and d2. Therefore,
% the definition of Fx,Fy here is k^2 times those in the single-frequency
clear;
clc;
close all
fclose('all');

addpath 'E:\My Drive\research\FDTD\Angora\ref\MatlabAnalysis'
% addpath '..\Tools'
addpath(genpath('C:\Users\Patrick\Desktop\PWS-SigmaToD'))




% this is meant to ignore border pixels who sometimes have edge artifacts
% bd = 0; %15;
% bdPer = 100; % the inner x%
bdPx = inf; %20; %inf;

dx = 10;

outputfullpath = 'E:\My Drive\research\FDTD\chromSTEM\20230315_NAi_055\output_glassNAi55\imaging';
num_of_samples = 1;

% number of plane-waves for each sample
number_of_planewaves_per_sample = 20; % 20;  % 20 is 10 for NAi >0, and x2 for each co and cross pol

if(number_of_planewaves_per_sample > 2)
    weight = [5.6220521601690017,2.3405721036739710,3.6562966125401654,4.0439940311519321,3.3673931131201436,2.9259610746657166,2.3976014830789016,1.7910853224061023,1.1088692300843212,0.40659940062234004]...
        .*[1,14,14,14,14,14,14,14,14,14]*1e-2/pi;
    weight = [weight,weight];
elseif(number_of_planewaves_per_sample == 2)
    weight = [1,1]/2;
else
    weight = 1;
end


NAc = sind(64.16);%1; %.98; %.63; %.98; %[.53, .66, .79, 0.98]  % RADIANS


% numerical aperture (NA) half angle
% Should match the collection angle in the Angora simulation  %asin(NA);  %NA=0.6
ap_angle = asin(NAc); % however, this in radians, and in the config file it's in degrees

% if true, the far-field file (.ff) is used to calculate the image again
use_far_field_file = true;

% normalize by the incident field intensity?
normalize_by_incident = true; %false;

if (use_far_field_file)

    % Magnification of the microscope : M
    % Sine of the angle w.r.t the z-axis, therefore the direction cosines w.r.t the x and y-axes are divided by this amount.
    % This results in the magnification of the image, because of the inverse relation between angles and distances (similar to the Fourier relation).
    % As the rays get more parallel to the axis, the image gets bigger.
    M = 20;
    image_extension_factor = 16;

    % obs_plane = 'down';
    obs_plane = 'up';

    % coordinates of the center of the image (in grid cells)
    % the coordinates are either w.r.t. the NFFFT origin, or the grid origin
    % (determined by "plot_wrt_grid_origin" below: If set, everything is w.r.t the grid origin, otherwise, the NFFFT origin)
    image_coord_x = 0;
    image_coord_y = 0;
    image_coord_z = 1e-6;%2e-6;

    % refractive index of the object space (THIS SHOULD ALWAYS BE 1)
    object_space_refr_index = 1;

    % refractive index of the image space (normalized by the object-side refractive index, since this is how it always appears in the formulas)
    image_space_refr_index = 1/1.517;

    % is phase-contrast applied to the image?
    phase_contrast = false;

    % take the dark-field image (suppress the incident field)?
    dark_field = false;
end



global_sample_index = 0;

intensity_tot_pixels_array = [];
intensity_unsca_pixels_array = [];

%% Loop through the samples
for sampleindex = 1:num_of_samples
    for planewaveindex = 1:number_of_planewaves_per_sample  % loop through the plane waves (for Kohler-illumination)

        farfieldfilename = fullfile(outputfullpath, ['Image_' num2str(global_sample_index) '_0.hd5.ff']);

        if (exist(farfieldfilename,'file'))
            if (use_far_field_file)

                % invoke the m-file that does the core operations for calculating the refocused field components
                calculate_refocused_field;

                fx=squeeze(f(1,:,:,:));
                fy=squeeze(f(2,:,:,:));
                fz=squeeze(f(3,:,:,:));
                E_pw_x=squeeze(E_pw(1,:,:,:));
                E_pw_y=squeeze(E_pw(2,:,:,:));
                E_pw_z=squeeze(E_pw(3,:,:,:));

                k_dim = 1;  % k dimension is the 1st
                x_dim = 2;  % x dimension is the 2nd
                y_dim = 3;  % y dimension is the 3rd
                center_pixel_x = floor(size(fx,x_dim)/2)+1;
                center_pixel_y = floor(size(fx,y_dim)/2)+1;
                %                 pixels_x = center_pixel_x+[-floor(num_pixels_x/2):floor((num_pixels_x-1)/2)];
                %                 pixels_y = center_pixel_y+[-floor(num_pixels_y/2):floor((num_pixels_y-1)/2)];
                %
                %                 % NEW
                %                 pixels_x = floor(length(x_range)/2)-2:floor(length(x_range)/2)+3;
                %                 pixels_y = floor(length(y_range)/2)-2:floor(length(y_range)/2)+3;
                %
                %                 % ALL
                %                 pixels_x = 1:length(x_range);
                %                 pixels_y = 1:length(y_range);

                %                 pixels_x = bd+1:length(x_range)-bd;
                %                 pixels_y = bd+1:length(y_range)-bd;

                %                 pixels_x = round((50-bdPer/2)*length(x_range)*.01)+1 : round((50+bdPer/2)*length(x_range)*.01);
                %                 pixels_y = round((50-bdPer/2)*length(y_range)*.01)+1 : round((50+bdPer/2)*length(y_range)*.01);

                if(isinf(bdPx))
                    pixels_x = 1:length(x_range);
                    pixels_y = 1:length(y_range);
                else
                    pixels_x = ceil(length(x_range)/2-bdPx/2):floor(length(x_range)/2+bdPx/2);
                    pixels_y = ceil(length(y_range)/2-bdPx/2):floor(length(y_range)/2+bdPx/2);
                end

                % Limit the image arrays to the pixels we are
                % interested in
                % fx,fy,fz: field "scattered" from inhomogeneities
                % E_pw_x,E_pw_y,E_pw_z: "unscattered" field, meaning
                % the Fresnel reflections from the layering
                fx_pixels = fx(:,pixels_x,pixels_y);
                fy_pixels = fy(:,pixels_x,pixels_y);
                fz_pixels = fz(:,pixels_x,pixels_y);
                E_pw_x_pixels = E_pw_x(:,pixels_x,pixels_y);
                E_pw_y_pixels = E_pw_y(:,pixels_x,pixels_y);
                E_pw_z_pixels = E_pw_z(:,pixels_x,pixels_y);

                if (~dark_field)
                    if (phase_contrast)
                        fx_pixels = fx_pixels + (1i)*E_pw_x_pixels; %#ok<UNRCH>
                        fy_pixels = fy_pixels + (1i)*E_pw_y_pixels;
                        fz_pixels = fz_pixels + (1i)*E_pw_z_pixels;
                    else
                        % Here, fx,fy,fz become the "total" field =
                        % "scattered"+"unscattered"
                        fx_pixels = fx_pixels + E_pw_x_pixels;
                        fy_pixels = fy_pixels + E_pw_y_pixels;
                        fz_pixels = fz_pixels + E_pw_z_pixels;
                    end
                end
            else
                hdf5_filename = fullfile(outputfullpath,['Image_' num2str(global_sample_index) '_0.hd5']); %#ok<UNRCH>

                kk = h5read(hdf5_filename,'/k_range');
                x_range = h5read(hdf5_filename,'/x_range');
                y_range = h5read(hdf5_filename,'/y_range');

                intensity_tot = h5read(hdf5_filename,'/intensity_tot');intensity_tot=permute(intensity_tot,[3,1,2]);
                intensity_unsca = h5read(hdf5_filename,'/intensity_unsca');intensity_unsca=permute(intensity_unsca,[3,1,2]);

                %                 pixels_x = floor(length(x_range)/2)-1:floor(length(x_range)/2)+3;
                %                 pixels_y = floor(length(y_range)/2)-1:floor(length(y_range)/2)+3;
                %                 pixels_x = bd+1:length(x_range)-bd;
                %                 pixels_y = bd+1:length(y_range)-bd;
                pixels_x = round((50-bdPer/2)*length(x_range)*.01)+1 : round((50+bdPer/2)*length(x_range)*.01);
                pixels_y = round((50-bdPer/2)*length(y_range)*.01)+1 : round((50+bdPer/2)*length(y_range)*.01);

            end

            if (planewaveindex==1)
                intensity_tot_pixels = zeros(length(kk),length(pixels_x),length(pixels_y));
                intensity_unsca_pixels = zeros(size(intensity_tot_pixels));
            end

            if (use_far_field_file)
                intensity_tot_pixels = intensity_tot_pixels + weight(planewaveindex)*(abs(fx_pixels).^2+abs(fy_pixels).^2+abs(fz_pixels).^2);
                intensity_unsca_pixels = intensity_unsca_pixels + weight(planewaveindex)*(abs(E_pw_x_pixels).^2+abs(E_pw_y_pixels).^2+abs(E_pw_z_pixels).^2);
            else
                intensity_tot_pixels = intensity_tot_pixels + weight(planewaveindex)*intensity_tot(:,pixels_x,pixels_y); %#ok<UNRCH>
                intensity_unsca_pixels = intensity_unsca_pixels + weight(planewaveindex)*intensity_unsca(:,pixels_x,pixels_y);
            end
        end

        %increment global sample index
        global_sample_index = global_sample_index + 1;
    end

    if (normalize_by_incident)
        intensity_tot_pixels = intensity_tot_pixels./intensity_unsca_pixels;
    end

    % Here, intensity_tot_pixels and intensity_unsca_pixels are 3D
    % arrays (1st:k, 2nd:x, 3rd:y)
    % Write the current image into the image array
    % (each column --1st. dim-- is a center spectrum)
    sample_dim = 4;
    intensity_tot_pixels_array = cat(sample_dim,intensity_tot_pixels_array,intensity_tot_pixels);
    intensity_unsca_pixels_array = cat(sample_dim,intensity_unsca_pixels_array,intensity_unsca_pixels);
end
%% Vis
img     = squeeze(mean(intensity_tot_pixels,1)); 
img_RGB = spectra2RGB(permute(intensity_tot_pixels,[2,3,1]),lambda*1e9);

img     = imrotate(img,90); 
img_RGB = imrotate(img_RGB,90); 
figure(11); 
subplot(121); imagesc(x_range*1e6,y_range*1e6,img);      colormap gray; axis equal; axis tight; xlabel('x (um)'); ylabel('y (um)'); title('Bright Field');  
subplot(122); imagesc(x_range*1e6,y_range*1e6,img_RGB);  colormap gray; axis equal; axis tight; xlabel('x (um)'); ylabel('y (um)');title('Bright Field -> CIE color matching');

% %%
% % figure;
% for ii = 1:num_of_samples
% % figure(ii); clf;
% %     subplot(121); imagesc(squeeze(mean(intensity_tot_pixels_array(:,:,:,ii),1))); axis equal; axis tight;  title('BF Image')
% %     subplot(122); imagesc(squeeze(std(intensity_tot_pixels_array(:,:,:,ii),[],1)));axis equal; axis tight;  title('\Sigma')
% %     fixfig
%
%     sigma_mean(ii) = median(makevec(std(intensity_tot_pixels_array(:,:,:,ii),[],1)));
%     D_all(ii,:) = SigmaToD(makevec(std(intensity_tot_pixels_array(:,:,:,ii),[],1)),1, .55, 0);
%     D_mean1(ii) = median(D_all(ii,:));
%     D_std1(ii) = std(D_all(ii,:));
%     sigma_std(ii) = std(makevec(std(intensity_tot_pixels_array(:,:,:,ii),[],1)));
% %     subplot(2,5,ii); imagesc(squeeze(std(intensity_tot_pixels_array(:,:,:,ii),[],1)),[0,.3]);axis equal; axis tight;  title('\Sigma'); fixfig
% end
% D_mean = SigmaToD(sigma_mean,1, .55, 0)
% D_std = D_mean-SigmaToD(sigma_mean-sigma_std,1, .55, 0)
%% Sigma2D_Patrick
load("E:\My Drive\research\FDTD\chromSTEM\ChromSTEMfile\ChromSTEM\RIproperties.mat")

cellRI = S2D.RIDefinition(ri_media,ri_chromatin,sigma_n);
sys = S2D.SystemConfiguration(cellRI,0.55,NAc*1.517,centerWavelength,true,true);
for ii = 1:num_of_samples

    data = squeeze(intensity_tot_pixels_array(:,:,:,ii));
    figure
    imagesc(squeeze(mean(data,1)))
    title('BF')
    figure;
    imagesc(squeeze(std(data,0,1)))
    colorbar
    title('\Sigma')


    [dOut, dCorrected, Nf_expected, lmax_corrected] = SigmaToD_AllInputs(squeeze(std(intensity_tot_pixels_array,1,1)), sys, 5e5, 2);
    figure
    imagesc(dCorrected,[1.8 3]);
    title('D');
    colorbar
    colormap jet

    figure
    histogram(dCorrected);
    title('D distribution')
end

% figure;
% %     subplot(221); imagesc(squeeze(mean(dn3b,3)),[1.33,1.37])
% %     subplot(221);imagesc(squeeze(mean(intensity_tot_pixels_array(:,:,:,ii)-1,1))',[-.05,.1]); axis equal; axis tight; colorbar; title('mean');
%     subplot(221);imagesc(squeeze(std(intensity_tot_pixels_array(:,:,:,ii)-1,[],1))',[0,.4]); axis equal; axis tight; colorbar; title('\Sigma');
% %     subplot(223); plot(squeeze(median(median(intensity_tot_pixels-1,2),3))); ylim([-.2,.2])
%     subplot(222);imagesc(reshape(D_all(ii,:),[size(intensity_tot_pixels_array,2),size(intensity_tot_pixels_array,3)])',[1.8,3.0]); axis equal; axis tight; colorbar; title('D')
%     subplot(223);histogram(D_all(D_all(:)>1.8),10); xlim([1.8,3.0]);
%
% %     title(sprintf('D = %4.2f/%4.2f',D_mean1(ii),D_mean(ii)));
%     set(gcf,'position',[809   295   956   683])
%             fixfig
%
% % return;
% %% Analyze D from polymer to compare to simulated
%
% D_acf = [];%[2.8248    2.7471    2.7268    2.7204    2.6860    2.6801    2.6661    2.6714    2.6805    2.6743];
% % %{
%
% fitL = round(150/(dx));
% for ai = 0 %:9
% %     progressbar(ai/9)
%     dn3b = GeomToMat(['G:\My Drive\research\wholeNucleusFDTD\20200221_fullspectra\input\random_medium_' num2str(ai) '.geom'],false);
%     [~,r,bnr] = acf(dn3b-mean(dn3b(:)),dx*1e-3);  %
%     figure(2); clf;
%         plot(log(r),log(bnr),'-'); hold on;
%         plot(log(r(2:fitL)),log(bnr(2:fitL)),'.'); fixfig
%     p = polyfit(log(r(2:fitL)),log(bnr(2:fitL)),1);
%     D_acf(ai+1) = p(1)+3;
% end
% %}
% %%
% fitL = round(150/(dx));
% p = polyfit(log(r(2:fitL)),log(bnr(2:fitL)),1);
% D_acf = p(1)+3;
% xPos = [0:519]*dx*1e-3; xPos = xPos-mean(xPos);
%  p = polyfit(log(r(2:fitL)),log(bnr(2:fitL)),1);
% figure(1); clf;
% %     subplot(221); imagesc(squeeze(mean(dn3b,3)),[1.33,1.37])
% %     subplot(221);imagesc(squeeze(mean(intensity_tot_pixels_array(:,:,:,ii)-1,1))',[-.05,.1]); axis equal; axis tight; colorbar; title('mean');
%     subplot(221);imagesc(xPos,xPos,squeeze(std(intensity_tot_pixels_array(:,:,:,ii)-1,[],1))',[0,.4]); axis equal; axis tight; colorbar; title('\Sigma');
%                 xlabel('x (um)');ylabel('y (um)');
% %     subplot(223); plot(squeeze(median(median(intensity_tot_pixels-1,2),3))); ylim([-.2,.2])
%     subplot(222);imagesc(xPos,xPos,reshape(D_all(ii,:),[size(intensity_tot_pixels_array,2),size(intensity_tot_pixels_array,3)])',[1.8,3.0]); axis equal; axis tight; colorbar; title('D')
%                 xlabel('x (um)');ylabel('y (um)');
%     subplot(224);hist(D_all(D_all(:)>1.8),100); xlim([1.8,3.0]); xlabel('D'); hold on ;plot(D_acf*ones(2,1),[0,15],'r--'); ylim([0,15]); title('D distribution'); ylabel('Count')
%     subplot(223); loglog(r,bnr,'-'); hold on;
%                   loglog(r(2:end),exp( polyval(p,log(r(2:end)))),'r:')
%                   loglog(r([2,fitL]),bnr([2,fitL]),'r.');
%
%                   ylim([.012,.679]*1e-4)
%                   xlabel ('r (um)');
%                   title('ACF of mass density B(r)');
%                   legend('ACF','Measured D','location','sw')
%
% %     title(sprintf('D = %4.2f/%4.2f',D_mean1(ii),D_mean(ii)));
%     set(gcf,'position',[809   295   956   683])
%     fixfigtimes
%
% %%
% warning off;
% D_mean = SigmaToD(sigma_mean,1, .55, 0)
% D_std = D_mean-SigmaToD(sigma_mean-sigma_std,1, .55, 0);

