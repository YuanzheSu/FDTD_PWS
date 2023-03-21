## FDTD_PWS
### Problems
#### 1. Bright field image of synthesized bright field image looks crapy
![](https://github.com/YuanzheSu/FDTD_PWS/blob/main/asset/brightfield.bmp)
<p align = "center">
Fig.1 - Bright field,M=20,OS=16
</p>

![](https://github.com/YuanzheSu/FDTD_PWS/blob/main/asset/fig1_M=1.bmp)
<p align = "center">
Fig.1x - Bright field,M=1,OS=1
</p>


#### 2. Sigma is 10x larger than predicted
Media is custom generated, and the calculated $\Sigma$ is way larger than theory value.
![](https://github.com/YuanzheSu/FDTD_PWS/blob/main/asset/sigma.bmp) 
<p align = "center">
Fig.2 - Spectra,M=20,OS=16
</p>

![](https://github.com/YuanzheSu/FDTD_PWS/blob/main/asset/fig2_M=1.bmp)
<p align = "center">
Fig.2x - Spectra,M=1,OS=1. Spectra dominated by green as can be seen from bright field.
</p>

Predicted $\Sigma$ is <0.1. 
I am not sure does problem 1 caused 2 or are they independent.

### Details
#### 1. Random Material Generation
See ```
code\FDTD_write_3D_random_material_Pat.m```


RM with $n = 1.38$, $\sigma_n=0.05$, $D=2.7$, $r_{min}=1$, $r_{max} = 100$
#### 2. FDTD
See ```FDTD\chromatin_on_glass_NAi55_a20.cfg``` Slightly modified from one of your previous config file

##### 2.1 Grid setup
$dx = 10 nm$, $X=Y=3 \mu m$, $Z = 4 \mu m$, 1 $\mu m$ thick glass on top, 3 $\mu m$ thick media on bottom.
![](https://github.com/YuanzheSu/FDTD_PWS/blob/main/asset/geom.png) 
Time step = 8000. Sufficient, see field value recording at center.
![](https://github.com/YuanzheSu/FDTD_PWS/blob/main/asset/fieldvaluerecording.bmp)
#### 2.2 Waveforms
500-700 nm.
#### 2.3 Incident
Used your method of using 20 plane waves to describe NAi=0.55
![](https://github.com/YuanzheSu/FDTD_PWS/blob/main/asset/incident.png)
#### 2.3 Recording
```
OpticalImages:
(
	{
	  output_data = ["intensity_tot","intensity_sca","intensity_unsca"];
	  num_of_lambdas = 30;
	  lambda_min = 500e-9;
	  lambda_max = 700e-9;
	  lambda_spacing_type = "k-linear";
	  ap_half_angle = 64.16; //53.14;
	  magnification = 1;
	  image_origin_z_in_cells = 100;  // glass-media interface
	  coll_half_space = "upper";
	}
);
```
#### 3. From .ff file to data cube
See ```code\fftosigma``` Slightly modified from your ```generatermsfrompolymer.m``` file.
The variables settings are below
```matlab
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

```
