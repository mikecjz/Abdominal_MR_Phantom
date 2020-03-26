function mixsamp = voximg2ksp(imPall,cmap,nval,opts)
% Convert phantom images to k-space
%   
%   INPUTS:
%       imPall [nx, ny, nz, 2] -> 4D phantom images generated by model2voximg.m
%       cmap   [nx, ny, nc, nz] -> coil sensitivity maps generated by gencmap.m
%       nval   [1 x 1]          -> noise level, generated from calcnoiselvl.m
%       opts   [struct]         -> structure for gridding non-Cartesian k-space with the NUFFT
%
%   OUTPUT:
%       mixsamp [nx, ny, nc, nz, nt] -> 4D phantom in k-space with phase accumulation and noise
%
% -----------------------------------------------------------------------------------------
% Realistic 4D abdominal phantom for magnetic resonance imaging
% Wei-Ching Lo
% wxl317@case.edu
% Case Western Reserve University
% April 2018
% -----------------------------------------------------------------------------------------

[nr,np] = size(opts.kx);
[~,~,npar,wfcomp,nframe] = size(imPall);
[~,~,nc,~] = size(cmap);
im2kspmask = true(nr,np);
    
mixsamp = zeros(nr,np,nc,npar,'single');

tstart = tic;
imW = permute(repmat(imPall(:,:,:,1),[1 1 1 nc]),[1 2 4 3]).*cmap;
imF = permute(repmat(imPall(:,:,:,2),[1 1 1 nc]),[1 2 4 3]).*cmap;

kspW = fft1n(imW,4);
kspF = fft1n(imF,4);

for c = 1:nc
    if strcmp(opts.trajectory,'cartesian')
        ksp2DW = squeeze(fft2n(kspW(:,:,c,:),1,2)); % one slice/partition
        ksp2DF = squeeze(fft2n(kspF(:,:,c,:),1,2));
    else
        for ipar = 1:npar
            ksp2DW(:,:,ipar) = embed(opts.G*kspW(:,:,c,ipar),im2kspmask);
            ksp2DF(:,:,ipar) = embed(opts.G*kspF(:,:,c,ipar),im2kspmask);
        end
    end
    
    % Apply phase accumulation
    mixksp = ksp2DW+ksp2DF.*repmat(exp(1i*2*pi*opts.FWshift*30*10^-6*(0:nr-1))',[1 np npar]);
    
    % Add noise
    mixsamp(:,:,c,:) = mixksp+nval*(randn(nr,np,npar)+1i*randn(nr,np,npar));
end

disp(['Simulate current frame in ' num2str(toc(tstart)) ' sec']);