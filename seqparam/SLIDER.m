
% -----------------------------------------------------------------------------------------
% Realistic 4D abdominal phantom for magnetic resonance imaging
% Wei-Ching Lo
% wxl317@case.edu
% Case Western Reserve University
% April 2018
% 
% 
% Created for SLIDER Sequence
% Junzhou Chen
% Cedars-Sinai Medical Center Biomedical Imaging Research Institute
% Junzhou.Chen@cshs.org
% April 2020

% -----------------------------------------------------------------------------------------

% Sequence Parameters: 'SLIDER' => Real-time multicontrast display for
% MR/RT

defseq.nPulseInGroupTR= 400; %number of pulses in groupTR
defseq.groupTR = 3200/1000; % sec
defseq.baseTR = 6/1000; % sec
defseq.baseTE = 2.12/1000; % sec
defseq.satTime = 100/1000; % sec Wait time after saturation recovery pulse
defseq.T2PrepTime = 66/1000; % sec Time between 2 90 pulses in T2prep module
defseq.gapTime = 700/1000; % sec gap time in the SLIDERsequence

defseq.isMultiFA = true;
defseq.FA = [5,10]; % degrees
defseq.FAStartLoc = [1/defseq.nPulseInGroupTR,326/defseq.nPulseInGroupTR];

defseq.isMultiPrep = true;
defseq.specialPrepLoc = [1/defseq.nPulseInGroupTR,301/defseq.nPulseInGroupTR,326/defseq.nPulseInGroupTR];


defseq.TRinprep = np*npar*nset; % # of TRs in one preparation
defseq.TRinex = 1; % # of TRs in one excitation
defseq.nprep = 1; % # of preparations(was 1)
defseq.phaseang = 0; % degrees
defseq.preppausetime = 0; % sec
defseq.dfrange = 0; % off-resonance in Hz
defseq.fatsat = 0; % fat saturation

% Preparation:
% 0: No prepatation
% 1: Inversion recovery
% 2: Spin echo
% 3: Saturation Recovery
% 4: T2prep (2 90 pulses in the same direction)
% 5: Gap
defseq.prepind = 0; % deafault no prep

defseq.specialPrep = [3,5,4]; % Special preperations for the SLIDER sequence. Use together with SpecialPrepLoc

% Define sampling pattern
% 'projinpar': projection in partition
% 'parinproj': partition in projection
% 'userdefined': user defined sampling pattern
defseq.samptype = 'projinpar';
if strcmp(sampmode,'demo')
    defseq.demosig = [1,100,200,300,315,350]; % selected contrast for phantom image
end