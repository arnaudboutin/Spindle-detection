function o_SS = ld_addSleepStage2spindles( i_SS, i_Infos , i_scoringFile, i_verb)
% 
% function SS_filtered = filter_spindles( i_SS, i_Infos )
% 
% i_SS: Spindle from mensen scripts
% i_Info: Information about spindles detection algorithm, electrodes,
%         acquisition
% 
% o_SS_filtered: Spindles with scoring label and IndChan label to know
% where spindle starts from.
% 
% Arnaud Bore: 20 Septembre 2015
%       Creation
% 
% Arnaud Bore: 26 Septembre 2015
%       Rename filter_spindles to addSleepStage2spindles
% 
if nargin < 4 
    i_verb = 1;
end

% Load marker file
if nargin < 3
    [FileName,PathName] = uigetfile('*.mat','Select marker file');
    markers = load([PathName,FileName]);
else
    markers = load(i_scoringFile);
end
% Acquisition time EEG (seconds)
eegLength = i_Infos.Recording.dataDim(2)/i_Infos.Recording.sRate;

% Scoring time (seconds)
scoringLength = length(markers.hypno_sieste{1,1})*markers.hypno_sieste{2,1};

if scoringLength < eegLength
    disp(['Scoring Length : ' num2str(scoringLength)]);
    disp(['EEG Length : ' num2str(eegLength)]);
    disp('Not enough scoring')
end
    
scoring = markers.hypno_sieste{1,1};

for nSp=1:length(i_SS)  % Loop spindles
    if i_verb
        disp(['Spindle: ' num2str(nSp)]);
    end
    for nRef=1:length(i_SS(1).Ref_Start) % Number of references
        % Get corresponding scoring Start and End spindle
        spStartSec = i_SS(nSp).Ref_Start(nRef) / i_Infos.Recording.sRate;
        spStartScoring = floor(spStartSec/markers.hypno_sieste{2,1})+1;

        spEndSec = i_SS(nSp).Ref_End(nRef) / i_Infos.Recording.sRate;
        spEndScoring = floor(spEndSec/markers.hypno_sieste{2,1})+1;

        if spStartSec==0
            i_SS(nSp).scoring(nRef) = nan;
        else
            if i_verb
                disp(['Ref: ' num2str(nRef)])
                disp(['spStartScoring: ' num2str(spStartScoring)])
            end
            if spStartScoring > size(scoring,2)
                if i_verb
                    disp('Spindle out of sleep scoring stages');
                end
                
                i_SS(nSp).scoring(nRef) = 0;
            else
                i_SS(nSp).scoring(nRef) = scoring(1,spStartScoring);
            end
        end
    end
end
o_SS = i_SS;

end

