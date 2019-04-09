function [] = plot_states_kappa(mat_All_mJointDist,nSubj,figure_directory)
    global m
    global tstart
    global tend
    global stepsize
    
    TMAX = size(mat_All_mJointDist{1},2);
    h = figure;
    hold on ,
    mKappaPlot = zeros(m, TMAX);
    for timeCounter  = 1:length(tstart:stepsize:tend)
        mJointDist = mat_All_mJointDist{timeCounter};
        for trial =1:1:TMAX
            ix = 1;
            for i = 1:m
                for j=1:m
                    mKappaPlot(j,trial) = mKappaPlot(j,trial) + mJointDist(ix,trial);
                    ix = ix+1;
                end
            end
        end
       subplot(2,4,timeCounter), imagesc(mKappaPlot); colormap hot
    end 
     figure_name = strcat('kappa_subject_',num2str(nSubj));
     savepath = strcat(figure_directory,figure_name,'.fig');
     saveas(h, savepath);
     close all;
end