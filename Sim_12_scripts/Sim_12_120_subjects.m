nSubj=120;
nRlz=30;

cd ../
addpath('/SPM8_PATH/')
addpath(genpath('../.'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_12(nSubj,['Sim_12_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
