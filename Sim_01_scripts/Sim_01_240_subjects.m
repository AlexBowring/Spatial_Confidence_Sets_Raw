nSubj=240;
nRlz=30;

cd ../
addpath('/SPM8_PATH/')
addpath(genpath('../.'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_01(nSubj,['Sim_01_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
