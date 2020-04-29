%% main_phantom_generation
%
% Preparation for Abdominal_MR_Phantom generation:
% (1) user-defined file name
% (2) sequence type
% (3) sampling trajectory
% (4) FOV and spatial resolution
% (5) sampling lines
% (6) respratory motion file and temporal resolution
% (7) coil sensitivity file and coil numbers
% (8) Fat-water chemical shift frequency
%
%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define tissue properties and respiratory motion.
% Generate voxelized volumetric mask with indexes.
% ------------------------------------------------
%
%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define sequence parameters.
% Generate signal evolution using Bloch simulation.
% Generate k-space data for each time frame.
% ------------------------------------------------
%
%%%%%%%%%% PART 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data reconstruction and gridding
% ------------------------------------------------
%
% -----------------------------------------------------------------------------------------
% Realistic 4D abdominal phantom for magnetic resonance imaging
% Wei-Ching (Tina) Lo
% wxl317@case.edu
% Case Western Reserve University
% July 2018
% -----------------------------------------------------------------------------------------
 
clear; close all;

addpath('/home/junzhou/Desktop/selectSlice_GUI')

% Define file name
savename = 'voxim15';

% Define sequence parameters:
% 'SpoiledGradientEcho' => T1-weighted image
% 'SpoiledGradientEchoWithFatSat' => Perfusion-weighted image
% 'InversionRecoveryLookLocker' => T1 mapping
% 'SingleSpinEcho' => T2 mapping
% 'SingleSpinEchoWithFatSat' => Diffusion-weighted image
% 'MultiEchoSpoiledGradientEcho' => Proton density fat fraction
sigtype = 'SLIDER';

% Define sampling trajectory:
% 'cartesian'
% 'radial'
% 'spiral'
samptraj = 'radial';

% Define 3D resolution and FOV
fov = 550; % field-of-view (mm)
mtx = 320; % matrix size
npar = 40; % # of partitions
slthick = 6; % slice thickness (mm)
nset = 1; % # of sets

% Define number of sampling lines
np_cartesian = 160; % # of Cartesian lines
np_radial = 100; % # of projections/spokes
np_spiral = 48; % # of spiral arms
sampmode = 'simple'; % sampling mode: 'demo', 'simple', 'eachEcho' 
% Note that "sampmode" will affect the simulation time (mins to hours, depending on the number of Echoes).
% Use 'eachEcho' option only when Echo-by-Echo simulation is nessasary

% Define 3D respiratory motion curve and temporal resolution (eq. each TR or each volume)
respmotion = 'respmov_long.mat';
tempres = 60; % temporal resolution (ms)
tempdur = 4200; % duration of each respiratory cycle (ms)
mscale = 1; % motion scale default 1
SImov = 13*mscale; % largest superior-inferior (SI) excursion (mm)
APmov = 6.5*mscale; % largest anterior-posterior (AP) excursion (mm)
LRmov = 2*mscale; % largest left-right (LR) excursion (mm)

% Define coil sensitivity map
coilmap = 'origcmap.mat';
nc = 18; % # of coils

% Define fat-water chemical shift
FWshift = 220; % water and fat are separated by approximately 440Hz in a 3T static field

% Load defined tissue properties: T1, T2, ADC, PDFF
tissueprop = tissueproperty;
tframe = floor(tempdur/tempres);

% Load ref images
load('/home/junzhou/AbdomenPhantom/phanimg_SLIDER_320_40_60ms.mat')
load('/home/junzhou/AbdomenPhantom/phanimg_long_SLIDER_320_40_60ms.mat')
load('/home/junzhou/AbdomenPhantom/phanimg_short_SLIDER_320_40_60ms.mat')
% combine phanimgs parfor experiments
% phanimg_all = cat(4,phanimg,phanimg_long,phanimg_short);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datapath = fileparts(mfilename('fullpath'));
cd(datapath)
addpath(genpath(datapath));

% Load defined tissue properties: T1, T2, ADC, PDFF
tissueprop = tissueproperty;

% respiratory motion pattern
load(respmotion); % respmov
tframe = floor(tempdur/tempres);
linv = genresp(respmov,tframe,[SImov APmov LRmov]);

% Convert mesh models to voxels (This step will take mins to hours, depending on the number of time frames)
xpts = -round((fov/fov)*mtx/2)+1:round((fov/fov)*mtx/2);
ypts = -round((fov/fov)*mtx/2)+1:round((fov/fov)*mtx/2); % x matrix size = y matrix size
zpts = linspace(round(-npar*slthick/3+10),round(npar*slthick/3+10),npar);
phanimg = mesh2model(tissueprop,linv,xpts,ypts,zpts);

% Show phantom images
showimg(phanimg);title('Indexed phantom images: axial plane');

save([savename '_phanimg.mat'],'phanimg','-v7.3')
fprintf('Phantom generation done\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare FFT/NUFFT
if strcmp(samptraj,'cartesian')
    np = np_cartesian;
elseif strcmp(samptraj,'radial')
    np = np_radial;
elseif strcmp(samptraj,'spiral')
    np = np_spiral;
end


opts = prepareNUFFT(mtx,np,samptraj,'linear_sorted','fov',fov,'FWshift',FWshift);

% Load in-vivo/simulated coil sensitivity maps
load(coilmap);
cmap = gencmap([mtx mtx npar],nc,origcmap);

% Sequence parameters
[seqparam,defseq] = setseqparam(sigtype,[np npar nset],sampmode);

%Overider demosig for SLIDER
defseq.demosig = 1:1:60000;
demosig = defseq.demosig;
% Bloch simulation
sigevo = gensigevo(tissueprop,seqparam);

% Convert voxels to phantom images (This step will take mins to hours, depending on the number of time frames)
nt = length(defseq.demosig);
[nr,~] = size(opts.kx);
% mixsamp = zeros(nr,np,nc,npar,nt,'single'); % will run out of memory for nt = 60,000

mixsampLin = zeros(nr,nc,nt)+1i*ones(nr,nc,nt);
kr = defseq.demosig;

sigma = 5; % std
kz = round(randn(1,nt)*sigma+20);
kz(kz>40) = 40;
kz(kz<1) = 1;

%random breathing pattern
respcycles = 200; %make sure it takes more than 8 minutes
cyclescat = generatecycletype(respcycles); % cycle catagories

[tpcat,tpOrder] = cycles2time(cyclescat); %time points catagories and order in their catagories. each timepoint is tempres apart;
respCat = zeros(1,nt);
refOrder = zeros(1,nt);
sigevoOrderAll = zeros(1,nt);
for i = 1:nt
    respCat(i) = tpcat(round(TRNumToTime(demosig(i))/tempres));
    refOrder(i) = tpOrder(round(TRNumToTime(demosig(i))/tempres));
    sigevoOrderAll(i) = demosig(mod(i-1,400)+1);
end

nval = 0.0016;
bigloopStart = tic;
parfor itp = 1:nt
    %     imPall =
    %     model2voximg(phanimg(:,:,:,mod(defseq.demosig(itp)-1,tframe)+1),sigevo(defseq.demosig(itp),:,:));
    %     % Ground truth images original
    s = tic;
    cat = respCat(itp);
    order = refOrder(itp);
    sigevoOrder = sigevoOrderAll(itp);
    
    switch cat
        case 1
            newOrder = order
            imPall = model2voximg(phanimg(:,:,:,newOrder),sigevo(sigevoOrder,:,:));
        case 2
            newOrder = order+70;
            imPall = model2voximg(phanimg_long(:,:,:,order),sigevo(sigevoOrder,:,:));
        case 3
            newOrder = order+175;
            imPall = model2voximg(phanimg_short(:,:,:,order),sigevo(sigevoOrder,:,:));
        otherwise
            error('Invalid respCat. Valid respCat {1,2,3}')
    end
%     imPall = model2voximg(phanimg(:,:,:,mod(round(TRNumToTime(demosig(itp))/tempres)-1,tframe)+1),sigevo(demosig(mod(itp-1,400)+1),:,:)); % Ground truth image
    
%         imPall(:,:,:,:,itp) = model2voximg(phanimg(:,:,:,mod(round(TRNumToTime(defseq.demosig(itp))/80)-1,tframe)+1),...
%             sigevo(defseq.demosig(itp),:,:)); % Ground truth images
%         if itp == 1
%             nval = calcnoiselvl(imPall,cmap);
%         end
    disp(itp)
    %     mixsamp(:,:,:,:,itp) = voximg2ksp(imPall,cmap,nval,opts); % k-space + noise
    
    %     ksp3D = voximg2ksppar(imPall,cmap,nval,opts); % k-space + noise
    
    if mod(itp-1,10) ==0 %navigator line
        %SLIDER Simulation only
%         cmdstr = 'opts_oneLine = prepareNUFFT(mtx,2,samptraj,''goldenAngle_360_tracked'',''fov'',fov,''FWshift'',FWshift,''GAtrace'',1);';
%         evalc(cmdstr);
        opts_oneLine = prepareNUFFT(mtx,2,samptraj,'goldenAngle_360_tracked','fov',fov,'FWshift',FWshift,'GAtrace',1);
        ksp2D = voximg2ksppar(imPall,cmap,nval,opts_oneLine,20); % k-space + noise
        mixsampLin(:,:,itp) = squeeze(ksp2D(:,1,:));
        
        kz(itp) = 20;
        kr(itp) = 20;
        disp(['loop in ',num2str(toc(s)),' sec'])
    else
        
        %SLIDER Simulation only
%         cmdstr = 'opts_oneLine = prepareNUFFT(mtx,2,samptraj,''goldenAngle_360_tracked'',''fov'',fov,''FWshift'',FWshift,''GAtrace'',itp);';
%         evalc(cmdstr);
        opts_oneLine = prepareNUFFT(mtx,2,samptraj,'goldenAngle_360_tracked','fov',fov,'FWshift',FWshift,'GAtrace',itp);
       
        ksp2D = voximg2ksppar(imPall,cmap,nval,opts_oneLine,kz(itp)); % k-space + noise
        mixsampLin(:,:,itp) = squeeze(ksp2D(:,1,:));
        disp(['loop in ',num2str(toc(s)),' sec'])
    end
end
disp(['SLIDER sim in ',num2str(toc(bigloopStart)),' sec'])

% save([savename '_mixsamp.mat'],'mixsamp','-v7.3')

save('../offline_ksp.mat','mixsampLin','-v7.3')
fprintf('Data acquisition done\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert 4D phantom k-space to images
reconimg = ksp2img(mixsamp,opts,cmap);
fprintf('Data reconstruction done\n');

% Show phantom images
showimg(reconimg(:,:,round(npar/2),:));colormap(gray);title('Reconstructed phantom images: axial plane')
showimg(imrotate(squeeze(reconimg(round(mtx/2),:,:,:)),90));colormap(gray);title('Reconstructed phantom images: frontal plane')


