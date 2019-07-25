nSubj=240;
nRlz=5;

cd /storage/maullz/Contour_Inference_2018/
addpath('/storage/essicd/spm8/')
addpath(genpath('/storage/maullz/Contour_Inference_2018/'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_11(nSubj,['Sim_11_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
