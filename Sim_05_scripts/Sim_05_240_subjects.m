nSubj=240;
nRlz=5;

cd /storage/maullz/Contour_Inference_2018/
addpath('/storage/essicd/spm8/')
addpath(genpath('/storage/maullz/Contour_Inference_2018/'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_05(nSubj,['Sim_05_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
