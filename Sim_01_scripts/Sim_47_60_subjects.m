nSubj=60;
nRlz=30;

cd /storage/maullz/Contour_Inference_2018/
addpath('/storage/essicd/spm8/')
addpath(genpath('/storage/maullz/Contour_Inference_2018/'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_47(nSubj,['Sim_47_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
