%% Sleep spindle detection
% Arnaud Boutin (CIAMS, Université Paris-Saclay)
% Adrien Conessa (CIAMS, Université Paris-Saclay)

%%%%%%
% Create "raw" and "output" folders
% The raw folder includes eeg files (.dat or .eeg, .vmrk, .vhdr) and scoring files (modify the suffix of the files)
% EEG files: Binary format and multiplexed 
% Run "nl_PPMIAO_SS_detection.m" to create sleep spindle matrices
%%%%%%



clc
redo=1; % 1 if you want to erase previous sps files

eegFolder ='';
outputeegFolder = ''; %num2str(type_info(idx_info))

scoringFolder='';

eeglab
close(gcf)
clc

load('channelsInfos.mat');

chan2process={'Pz'}; 
pos_chan2pro=find(ismember([ChanInfos.labels],chan2process));


allFiles = dir([ eegFolder '*.eeg']);

% Spindle Detection %
% ----------------- %
for iFile=subject_to_process-4
    
    tic
    
    o_name = allFiles(iFile).name;
    eeg_name = [eegFolder allFiles(iFile).name];
    vhdr_name = [allFiles(iFile).name(1:end-3) 'vhdr'];
    
    out_name = ['o_' allFiles(iFile).name(1:end-3) 'mat'];
    
    fprintf(['\nSubject: ' num2str(iFile+4) '\n'])
    subject=num2str(iFile+4);

    if iFile+4<10
        subject_str=['00' subject];
    else
        subject_str=['0' subject];
    end
  
    
    % Reset channels info for each participant
    clear currentChanInfos
    

    if exist([outputeegFolder out_name],'file') && redo==0
        disp(['Already done : ' out_name]);
        load([outputeegFolder out_name])
    else
        
        disp(['Analysis: ' allFiles(iFile).name]);

        try % Read Scoring file
            i_scoringFile = [scoringFolder filesep 'hypno_S' subject_str '_nap.mat']; % scoring file suffix %  'score_test' 
            load (i_scoringFile)
        catch
            disp(['Error no scoring file found : ' allFiles(iFile).name])
        end

        sample_window=[]; % if you want to specify the eeg window you want to use, otherwise leave empty

        % load EEG data
        EEG = pop_loadbv(eegFolder, vhdr_name, sample_window, pos_chan2pro); % specify channels to load 
         
        for nCh=1:length(EEG.chanlocs)
            disp(['Channel : ' EEG.chanlocs(nCh).labels]);
        end

        % Info initialization
        Info.Recording.dataDim = size(EEG.data);
        Info.Recording.sRate = EEG.srate;

        % get the default settings for spindle detection
        try
            Info = swa_getInfoDefaults(Info, 'SS');
        catch
            disp(['Error swa_getInfoDefaults function: ' allFiles(iFile).name])
        end

        % Check here if trouble with multiple ch detection
        for nChan =length(EEG.chanlocs):-1:1  % number of channels
            Ind = find(strcmp([ChanInfos.labels],EEG.chanlocs(nChan).labels));
            if isempty(Ind)
                disp(['Error no ' EEG.chanlocs(nChan).labels ' ch found in channel infos'])
            else
                currentChanInfos(1,nChan) = ChanInfos(Ind);
            end
        end


        Info.Electrodes = currentChanInfos;
        if size(Info.Electrodes,2)>1
            [~,sortElec] = sort(convertStringsToChars([Info.Electrodes(:).labels]));
            Info.Electrodes = Info.Electrodes(sortElec);
        else
            sortElec=1;
        end
            

        Data.SSRef = EEG.data(sortElec,:);
            
        % Detection Spindle Mensen script
        [Data, Info, SS_Core] = ld_nl_swa_FindSSRef(Data, Info, hypno_sieste);

        % add sample idx before sleep time, removed previously
        if ~isempty (sample_window)
            for idx_sps=1:numel(SS_Core)

                idx_no_0 = find(SS_Core(idx_sps).Ref_Start~=0);
                add_sample= zeros(1,numel(SS_Core(idx_sps).Ref_Start));
                add_sample(idx_no_0)=sample_window(1)-1;


                SS_Core(idx_sps).Ref_Start = SS_Core(idx_sps).Ref_Start + add_sample;
                SS_Core(idx_sps).Ref_End = SS_Core(idx_sps).Ref_End + add_sample;

            end
        end

            
            try % Filter spindles depending on stage scoring
                SS = ld_addSleepStage2spindles(SS_Core, Info, i_scoringFile);
            catch
                disp(['Error ld_addSleepStage2spindles function: ' allFiles(iFile).name])                
            end
            
            try % Save output
                Size_file = whos('Size','SS'); % SS files need to be < 2 GB on Matlab v7 and sup
                    if Size_file.bytes < 2e+09
                        swa_saveOutput(Data, Info, SS, [outputeegFolder out_name], 0, 0);
                        disp(['Save ' out_name]);
                    else
                        disp('SS file bigger than 2GB');
                        swa_saveOutput_bigfile(Data, Info, SS, [outputeegFolder out_name], 0, 0);
                        disp(['Save ' out_name]);
                    end
            catch
                disp(['Error swa_saveOutput function: ' allFiles(iFile).name])                
            end
        
        clear Data EEG SS SS_Core Info_input Ind o_name i_marker name i_struct_marker currentChanInfos Info Size_file
    end
    
   

    clear Data Info SS SS_Core
    
    toc
end
