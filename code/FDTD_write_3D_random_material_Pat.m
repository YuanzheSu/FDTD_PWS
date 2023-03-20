clear;
clc;
close all
%%
% n0_array = 1.517;

% lc_array = [20,22.8,26.1,29.8,34,38.9,44.4,50.7,57.9,66.2,75.6,86.3,98.6,113,129,147,168,192,219,250]/20;
% lc_array = [22.8,29.8,38.9,50.7,66.2,86.3,113,147,192,250]/20;


% lc_array = 0/20;
mean_n  = 1.38;
sigma_n = 0.05;
D       = 2.7; %2.7;
r_min   = 1;
% r_max   = 251/20;
r_max   = 100;  % aya
% L       = 200;

%% extent of the random medium

% grid spacing in the fdtd grid (in meters)
dx = 10e-09;
x = 1.5e-6;
y = 1.5e-6;
z = 2e-6;

M = round(x/dx);     % x dimension 225 969
N = round(y/dx);     % y dimension 225 969
P = round(z/dx);    % z dimension -- should be greater than roughly (mean_thickness+6*sigma_h) 151 61

L = 2*P+1;

%%
num_of_samples = 1;


print_sample_number = true;
show_media_as_image = false;%true; %false;
append_medium_index = true;



% inputbasepath = ['C:\Users\Aya\Google Drive\SimulationCode\FDTD\RandomMedia\20171121_PWS_2piece_changingNAi\NAi0.55\D' num2str(D) '\input\'];%'./';
inputbasepath = 'E:\My Drive\research\FDTD\chromSTEM\20230215_calibration\';
directory = 'input_1\';

samplefilename = 'random_medium_';

if (~append_medium_index)
    warning('Random media indices are not appended to file names. This could result in overwriting.');
end
inputpath = [inputbasepath directory];


%%
% we use a global sample index to index all the samples, for ease of simulation in fdtd
global_sample_index = 0;



% disp([num2str(num_of_samples) ' samples of the random medium (rough surface) to be calculated for ' ...
%     'n0=' num2str(n0_array) ', lc=' num2str(lc_array) ]);

% % std. dev. of thickness fluctuation in grid cells
% sigma_h_cell = sigma_h/dx;
% % correlation length in grid cells
% corrlen = lc/dx;

for sample_index=1:num_of_samples
    %       Material  =  generate_corr_3D(sigma_n,lc,M,N,P,2)+mean_n; % material in middle
    %       bnr = ACF_fractal_2piece(sigma_n,D,r_min,r_max,L)';
    %       bnr = ACF_fractal_2piece_aya(sigma_n,D,r_min,r_max,L)';

    bnr = ACF_fractal_analytical(sigma_n,D,r_min,r_max,L)';

    %       bnr = ACF_fractal_2point(sigma_n,D,r_min,r_max,L)';

    Material = Bnr2medium(bnr,L)+mean_n;


    %       Material = Material(1:M, 1:N, 1:P);
    %       [~, ~, bnr_g, ~] = acf(Material-mean_n,1);
    %       figure;plot(bnr);hold on;plot(bnr_g,'r');
    %       xlabel('r(nm)');ylabel('Bnr')
    %     figure;loglog([0:L]+eps,bnr);hold on;loglog([0:L],bnr,'bx'); hold on;plot(r_min*ones(2,1),ylim,'k--');plot(r_max*ones(2,1),ylim,'k--')
    %     r = [0:L]; % for 2 point
    %     rFit = r>5*(r_min) & r<(r_max/5);
    %     p1 = polyfit(log(r(rFit)),log(bnr(rFit)'),1) ;
    %     p2 = polyfit(log(r(rFit)),real(log(bnr_g(rFit)')),1);
    %     [D, 3+p1(1), 3+p2(1)]

    Material = Material(1:M, 1:N, 1:P);

    if (append_medium_index)
        fullfilename = [inputpath samplefilename num2str(global_sample_index) '.geom'];
    else
        fullfilename = [inputpath samplefilename '.geom']; %#ok<*UNRCH> 
    end

    clear geomfile
    geomfile = fopen(fullfilename,'w');    %open file for writing

    % create the uniform slab

    eps_r = Material.^2;

    if (show_media_as_image)
        %         y=1:size(d_er,2);
        %         z=1:size(d_er,3);
        %         figure;imagesc(z*dx*1e9,y*dx*1e9,squeeze(sqrt(d_er(round(M/2),:,:))));axis image xy;colorbar;
        %         xlabel('z (\mum)');ylabel('y (\mum)');
        figure;
        im=imagesc(flipud(squeeze(eps_r(1,:,:)).'));
        hold on;
        for i=1:M
            set(im,'CData',flipud(squeeze(eps_r(i,:,:)).'));
            pause(.01);
        end
    end

    % write the dimensions first
    fwrite(geomfile,M,'int');
    fwrite(geomfile,N,'int');
    fwrite(geomfile,P,'int');
    % then write the array itself
    for k=1:size(eps_r,3)    % z-dimension written last
        fwrite(geomfile,eps_r(:,:,k),'double'); % write the x and y dimensions (in this order)
    end

    fclose('all');

    if (print_sample_number)
        disp(['Sample ' num2str(sample_index) ' written to ' fullfilename]);
    end


    %increment global sample index
    global_sample_index = global_sample_index+1;
end

