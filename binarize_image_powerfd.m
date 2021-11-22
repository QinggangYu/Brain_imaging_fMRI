load('/home/qinggang/research/MIDUS_refresher/resting_conn/powerfd_badframe_72.mat')
nrun = size(bad_f_clean, 2); % enter the number of runs here
jobfile = {'/home/qinggang/research/MIDUS_refresher/rsfmri/binarize_image_job.m'};
load('/home/qinggang/research/MIDUS_refresher/rsfmri/filenames_resting_powerfd.mat')
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    curr_file = clean_sub_fname(crun, :);
    %sub_id = curr_file(57:end);
    inputs{1, crun} = cellstr(fullfile(curr_file, 'wau_mean_powerfd.nii'));
    inputs{2, crun} = cellstr(curr_file);
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});