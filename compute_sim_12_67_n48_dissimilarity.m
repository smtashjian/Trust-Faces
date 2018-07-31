clear
clc
subids = {'101', '102', '103', '105', '106','107', '108', '109', '110', '111', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '125', '126', '127', '128', '129', '130', '131', '132', '134', '135', '137', '138', '140', '141', '142', '143','144','145', '146', '147', '148', '149', '150', '151', '152','153', '156'}; 
trials = {'01','02','03','04','05', '06', '07', '08', '09', '10', '11', '12','13', '14', '15', '16', '17', '18', '19', '20'};
%c = 1;
 simMat = zeros(numel(subids), 2);
 
 
for sub = 1:length(subids);
    data_path = '/u/project/agalvan/data/FSS_W2/Faces/Trust/RSA/';
    c = 1;
      dat_12_ffaneuro=[];
        dat_12_workingmem=[];
        dat_67_ffaneuro=[];
        dat_67_workingmem=[];
    for trl = 1:length(trials);       
        ds_12_ffaneuro = cosmo_fmri_dataset([data_path sprintf('FSS_%s_12.%s.feat/stats/zstat1.nii.gz', subids{sub}, trials{trl})], 'mask', [data_path sprintf('FSS_%s_12.%s.feat/reg/reg_amygdala_FSS_%s.nii.gz', subids{sub}, trials{trl}, subids{sub})]);
        ds_12_workingmem = cosmo_fmri_dataset([data_path sprintf('FSS_%s_12.%s.feat/stats/zstat1.nii.gz', subids{sub}, trials{trl})], 'mask', [data_path sprintf('FSS_%s_12.%s.feat/reg/reg_Juelich_centromedial_FSS_%s.nii.gz', subids{sub}, trials{trl}, subids{sub})]);
        ds_67_ffaneuro = cosmo_fmri_dataset([data_path sprintf('FSS_%s_67.%s.feat/stats/zstat1.nii.gz', subids{sub}, trials{trl})], 'mask', [data_path sprintf('FSS_%s_67.%s.feat/reg/reg_amygdala_FSS_%s.nii.gz', subids{sub}, trials{trl}, subids{sub})]);
        ds_67_workingmem = cosmo_fmri_dataset([data_path sprintf('FSS_%s_67.%s.feat/stats/zstat1.nii.gz', subids{sub}, trials{trl})], 'mask', [data_path sprintf('FSS_%s_67.%s.feat/reg/reg_Juelich_centromedial_FSS_%s.nii.gz', subids{sub}, trials{trl}, subids{sub})]);
        
        %dat_12_ffa = zeros(numel(ds_12_ffa.samples), numel(trials));
        %dat_67_ffa = zeros(numel(ds_67_ffa.samples), numel(trials));
        %dat_12_cingpregen = zeros(numel(ds_12_cingpregen.samples), numel(trials));
        %dat_67_cingpregen = zeros(numel(ds_67_cingpregen.samples), numel(trials));
        
        %dat_12_ffa=zeros(numel(ds_12_ffa.samples),numel(trials));
        %dat_12_cingpregen=zeros(numel(ds_12_cingpregen.samples),numel(trials));
        %dat_67_ffa=zeros(numel(ds_67_ffa.samples),numel(trials));
        %dat_67_cingpregen=zeros(numel(ds_67_cingpregen.samples),numel(trials));
        
        dat_12_ffaneuro(:,c) = ds_12_ffaneuro.samples';
        dat_12_workingmem(:,c) = ds_12_workingmem.samples';
        dat_67_ffaneuro(:,c) = ds_67_ffaneuro.samples';
        dat_67_workingmem(:,c) = ds_67_workingmem.samples';
        disp(sprintf('Finished with sub FSS%s, trial %s', subids{sub},trials{trl}));
        c = c+1;
    end
    dat_all_ffaneuro = [dat_12_ffaneuro, dat_67_ffaneuro];
    dat_all_ffaneuro = dat_all_ffaneuro(any(dat_all_ffaneuro,2), :);
    %dat_12_ffa = dat_12_ffa(any(dat_all_ffa,2),:);
    %dat_67_ffa = dat_67_ffa(any(dat_all_ffa,2),:);
    z_all_ffaneuro = zscore(dat_all_ffaneuro,0,1);
    ffaneuro_sim_mat = corr(z_all_ffaneuro, 'Type', 'Pearson');
    ffaneuro_12_67_sim = ffaneuro_sim_mat((1:20), (1:20));
    
    dat_all_workingmem = [dat_12_workingmem, dat_67_workingmem];
    dat_all_workingmem = dat_all_workingmem(any(dat_all_workingmem,2), :);
    z_all_workingmem = zscore(dat_all_workingmem,0,1);
    workingmem_sim_mat = corr(z_all_workingmem, 'Type', 'Pearson');
    workingmem_12_67_sim = workingmem_sim_mat((1:20), (1:20));
    
    
    %mean_12_ffa = mean(dat_12_ffa, 2);
    %mean_67_ffa = mean(dat_67_ffa, 2);
    
    %z_12_ffa = zscore(mean_12_ffa);
    %z_67_ffa = zscore(mean_67_ffa);
    
    %sim_ffa = corr(z_12_ffa,z_67_ffa);
    
    sim_f_ffaffaneuro = (0.5)*log((1+ffaneuro_12_67_sim)./(1-ffaneuro_12_67_sim));
    sim_f_ffaffaneuro = sim_f_ffaffaneuro(tril(sim_f_ffaffaneuro, -1) ~= 0); %sim_f_ffaffaneuro = sim_f_ffaffaneuro(:);
    sim_ffaneuro = mean(sim_f_ffaffaneuro);
    
    sim_f_workingmem = (0.5)*log((1+workingmem_12_67_sim)./(1-workingmem_12_67_sim));
    sim_f_workingmem = sim_f_workingmem(:);
    sim_workingmem = mean(sim_f_workingmem);
    
    simMat(sub, 1) = sim_ffaneuro;
    simMat(sub, 2) = sim_workingmem;
    
end    
