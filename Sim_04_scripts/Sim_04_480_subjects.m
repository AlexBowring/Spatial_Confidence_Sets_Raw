nSubj=480;
nRlz=30;

cd /storage/maullz/Contour_Inference_2018/
addpath('/storage/essicd/spm8/')
addpath(genpath('/storage/maullz/Contour_Inference_2018/'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_04(nSubj,['Sim_04_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
