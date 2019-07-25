nSubj=240;
nRlz=5;

cd ../
addpath('/SPM8_PATH/')
addpath(genpath('../.'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_05(nSubj,['Sim_05_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
