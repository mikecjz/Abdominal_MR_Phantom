function sigevo = gensigevo(tissueprop,seqparam)

% Generate signal evolution for each tissue with specified sequence parameters
%   
%   INPUTS:
%       tissueprop [struct] -> tissue property generated by tissueproperty.m
%       seqparam [struct] -> sequence parameter generated by demo_seqparam_release.m
%
%   OUTPUT:
%       sigevo [Nex x cnt x 2] -> signal evolution
%                                 Nex = # of excitation
%                                 cnt = # of tissues 
%                                 water(1) and fat(2) signal
%
% -----------------------------------------------------------------------------------------
% Realistic 4D abdominal phantom for magnetic resonance imaging
% Wei-Ching Lo
% wxl317@case.edu
% Case Western Reserve University
% April 2018
% -----------------------------------------------------------------------------------------

r = tissueprop;
FF = [tissueprop.PDFF];
ADC = [tissueprop.ADC];

FA = seqparam.FAmat;
TR = seqparam.TRmat;
TE = seqparam.TEmat;
prep = seqparam.prep;
ti = seqparam.ti;
phase = seqparam.phase;
df = seqparam.dfrange;
spoiler = seqparam.spoiler;
perfT1 = seqparam.perfT1;
Nex = seqparam.Nex;
Necho = seqparam.Necho;
cnt = numel(r);


if isfield(seqparam,'bval') && any(seqparam.bval > 0)
    bval = seqparam.bval;
end

if isfield(seqparam,'satTime')
    satTime = seqparam.satTime;
end

if isfield(seqparam,'T2PrepTime')
    T2PrepTime = seqparam.T2PrepTime;
end

if isfield(seqparam,'gapTime')
    gapTime = seqparam.gapTime;
end

%% Bloch simulation of 2D QALAS signal
% Spoiler gradient
Nspins = 100;
phiRange = linspace(0,2*pi,Nspins);
Z = zeros(3,3,Nspins);
for iph = 1:Nspins
    Z(:,:,iph) = zrot(phiRange(iph));
end

MX = zeros(cnt,Nex*Necho,'single'); % real valued, Mx
MY = zeros(cnt,Nex*Necho,'single'); % real valued, My
MZ = zeros(cnt,Nex*Necho,'single');

R180 = yrot(pi);
Rx90 = xrot(pi/2);

rfPulseDur = 0.8/1000;
tidefault = 20.64/1000;
for idf = 1:numel(seqparam.dfrange)
    fprintf('df = %d\n',seqparam.dfrange(idf));
    for ir = 1:cnt
        
        
%         if ir == 2
%             disp(ir)
%         elseif ir == 11
%             disp(ir)
%         end
        
        T1 = r(ir).T1/1000;
        T2 = r(ir).T2/1000;
        
        M0 = repmat([0;0;1],[1 Nspins]);
        M = zeros(3,Nspins,Nex,Necho);
        Mtmp = M0;        
        rf = 0; % initial rf phase
        
        for n = 1:Nex
            if any(perfT1(:)>0)
                if ir == 7
                    T1 = perfT1(floor((n-1)/size(prep,1))+1,4)/1000; % liver
                elseif ir == 11
                    T1 = perfT1(floor((n-1)/size(prep,1))+1,3)/1000; % vein
                elseif ir == 14
                    T1 = perfT1(floor((n-1)/size(prep,1))+1,2)/1000; % artery
                end
            end
            Mtmp(1:2,:) = 0;
            TErelax = TE(n,:);
            relaxTime = TR(n)-TErelax(1)-rfPulseDur/2;
            if TErelax < 0, error(['TErelax ' num2str(n) ' < 0']); end
            
%             if n == 325
%                 disp(n)
%             elseif n == 326
%                 disp(n)
%             elseif n == 400
%                 disp(n)
%             elseif n == 401
%                 disp(n)
%             end
            
            switch prep(n)
                
                case 0 % no prep
                    if relaxTime < 0, error(['relaxTime ' num2str(n) ' < 0']); end
                    Mtmp = freeprecspin(Mtmp,relaxTime,T1,T2,df,Nspins);
                case 1 % Inversion recovery
                    tiFill = ti(n) - tidefault;
                    if tiFill < 0, error(['tiFill ' num2str(n) ' < 0']); end
                    relaxTime = relaxTime - tidefault - tiFill;
                    if relaxTime < 0, error(['relaxTime ' num2str(n) ' < 0']); end
                    
                    Mtmp = freeprecspin(Mtmp,relaxTime,T1,T2,df,Nspins);
                    Mtmp = R180*Mtmp;
                    Mtmp = freeprecspin(Mtmp,tidefault,T1,T2,df,Nspins);
                    for iph = 1:Nspins
                        Mtmp(:,iph) = Z(:,:,iph)*Mtmp(:,iph);
                    end
                    if tiFill > 0
                        Mtmp = freeprecspin(Mtmp,tiFill,T1,T2,df,Nspins);
                    end
                case 2 % Spin echo
                    relaxTime = relaxTime - TErelax;%t2te(n);
                    if relaxTime < 0, error(['relaxTime ' num2str(n) ' < 0']); end
                    Mtmp = freeprecspin(Mtmp,relaxTime,T1,T2,df,Nspins);
                    Mtmp = Rx90*Mtmp;
                    Mtmp = freeprecspin(Mtmp,TErelax/2,T1,T2,df,Nspins);
                    TErelax = TErelax/2;
                case 3 % Saturation Recovry
                    Mtmp = freeprecspin(Mtmp,relaxTime,T1,T2,df,Nspins); %This step is perhaps unneccessary for SR? 
                    
                    Mtmp = Rx90*Mtmp;
                    Mtmp = freeprecspin(Mtmp,satTime,T1,T2,df,Nspins);
                    
                    Mtmp(1:2,:) = 0;% Assume spoiler
                case 4 % T2 prep (2 same direction 90 pulses)
                    Mtmp = freeprecspin(Mtmp,relaxTime,T1,T2,df,Nspins); %This step is perhaps unneccessary for T2prep? 
                    Mtmp = Rx90*Mtmp;
                    
                    Mtmp = freeprecspin(Mtmp,T2PrepTime,T1,T2,df,Nspins);
                    
                    Mtmp = -1*Rx90*Mtmp;
                    
                    Mtmp(1:2,:) = 0;% Assume spoiler
                case 5 % Gap
                     Mtmp = freeprecspin(Mtmp,gapTime,T1,T2,df,Nspins);
                otherwise
                    error('invalid prep pulse')
            end
            
            % Readout
            Mtmp = freeprecspin(Mtmp,rfPulseDur/2,T1,T2,df,Nspins); % RF
            Mtmp = yrot(FA(n))*Mtmp;
            for p = 1:Necho
                Mtmp = freeprecspin(Mtmp,TErelax(p),T1,T2,df,Nspins); % TE
                M(:,:,n,p) = Mtmp;
            end
            inc = phase(n);
            rf = mod(rf+inc,2*pi);
            Mtmp = zrot(rf)*Mtmp;
                
            if spoiler(n) == 1
                for iph = 1:Nspins
                    Mtmp(:,iph) = Z(:,:,iph)*Mtmp(:,iph);
                end
            end
            
        end
        % Average signal over all spins in voxel
        Mavg = squeeze(mean(M,2));
        MX(ir,:) = Mavg(1,:);
        MY(ir,:) = Mavg(2,:);
        MZ(ir,:) = Mavg(3,:);
    end % number of T1 and T2 combinations
    MX = permute(MX,[2 1]);
    MY = permute(MY,[2 1]);
    MZ = permute(MZ,[2 1]);
    if exist('bval','var')
        sigdata(:,(idf-1)*cnt+1:idf*cnt) = (MX + 1i*MY).*exp(-bval'*ADC);
    else
        sigdata(:,(idf-1)*cnt+1:idf*cnt) = (MX + 1i*MY);
    end
end

% Plot off-resonance at 0Hz
figure;hold on;
plot(abs(sigdata(:,2)),'b','linewidth',2); % fat
plot(abs(sigdata(:,11)),'g','linewidth',2); % vein
plot(abs(sigdata(:,7)),'k','linewidth',2); % liver
plot(MZ(:,1),'r-'); % MZ
legend('Fat','Vein','Liver','Mz')
drawnow;

% Generate signal for water and fat
if seqparam.fatsat == 1
    sigevo = zeros([Nex*Necho cnt],'single');
    sigevo(:,:,1) = reshape(sigdata.*repmat((1-FF),[Nex*Necho 1]),[Nex*Necho cnt]);
    sigevo(:,:,2) = zeros(Nex*Necho,cnt,'single');
else
    sigevo = zeros([Nex*Necho cnt 2],'single');
    sigevo(:,:,1) = reshape(sigdata.*repmat((1-FF),[Nex*Necho 1]),[Nex*Necho cnt]);
    sigevo(:,:,2) = reshape(sigdata.*repmat(FF,[Nex*Necho 1]),[Nex*Necho cnt]);
end

fprintf('Signal simulation done\n');