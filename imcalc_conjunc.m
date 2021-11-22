% List of open inputs
nrun = 264; % enter the number of runs here
jobfile = {'/home/qinggang/research/MIDUS_refresher/resting_conn/imcalc_conjunc_job.m'};
vol_f = '/home/qinggang/research/MIDUS_refresher/rsfmri/prod_bi50_powerfd.nii,1';
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    mask_f = sprintf('power4d_3mm.nii,%d', crun);
    mask = fullfile('/home/qinggang/research/MIDUS_refresher/resting_conn/', mask_f);
    inputs{1, crun} = {vol_f; mask};
    inputs{2, crun} = sprintf('junc_power%03d.nii', crun);
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
