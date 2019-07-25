nSubj=120;
nRlz=30;

cd ../
addpath('/SPM8_PATH/')
addpath(genpath('../.'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_02(nSubj,['Sim_02_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
