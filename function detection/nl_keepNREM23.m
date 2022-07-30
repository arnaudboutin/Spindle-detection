function Data = nl_keepNREM23(Data, hypno_sieste)

% Replace non-NREM23 sample by NaN value

NREM23_samples = zeros(1,numel(Data));


NREM2_samples = floor([hypno_sieste{3}.NREM2.onset*250 ; (hypno_sieste{3}.NREM2.onset + hypno_sieste{3}.NREM2.duration) *250]);
NREM3_samples = floor([hypno_sieste{3}.NREM3.onset*250 ; (hypno_sieste{3}.NREM3.onset + hypno_sieste{3}.NREM3.duration) *250]);
bad_samples = floor([hypno_sieste{3}.BadIntervals.onset*250 ; (hypno_sieste{3}.BadIntervals.onset + hypno_sieste{3}.BadIntervals.duration) *250]);
bad_samples(bad_samples==0)=1;  % if the bad intervals begin at t=0


for NREM2_wd=1:size(NREM2_samples,2)
    
    NREM23_samples(NREM2_samples(1,NREM2_wd):NREM2_samples(2,NREM2_wd))=1;
    
end

for NREM3_wd=1:size(NREM3_samples,2)
    
    NREM23_samples(NREM3_samples(1,NREM3_wd):NREM3_samples(2,NREM3_wd))=1;
    
end

for bad_wd=1:size(bad_samples,2)
   
        
    NREM23_samples(bad_samples(1,bad_wd):bad_samples(2,bad_wd))=0;
    
end



if sum(NREM23_samples)~=0
    Data(NREM23_samples==0) = nan;
end
