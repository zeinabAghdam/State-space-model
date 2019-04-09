function [mAvgKappa] = computeKappaLikelihood( mat_all_subjects_mJointDist, subj, numDecible, figure_directory_2, ...
    numConfiguration,TMAX)
% compute the kappa and mean matrix.
global tstart
global stepsize
global tend
global m
mAvgKappa = zeros(m,TMAX);
mAvgMu    = zeros(m,TMAX);

    mJointDist = mat_all_subjects_mJointDist{1};
    
    h = figure; hold on ,
    for timeCounter  = 1:length(tstart:stepsize:tend)
        mKappaPlot = zeros(m, TMAX);
        mMeanPlot = zeros(m,TMAX);        
        for trial =1:1:TMAX
            ix = 1;
            for i = 1:m
                for j=1:m
                    mKappaPlot(j,trial) = mKappaPlot(j,trial) + mJointDist(ix,trial);
                    mMeanPlot(i,trial) = mMeanPlot(i,trial)+ mJointDist(ix,trial); 
                    ix = ix+1;
                end
            end
        end
        %subplot(3,4,timeCounter), imagesc(mMeanPlot); colormap hot
        imagesc(mKappaPlot); colormap hot
        mAvgKappa = mAvgKappa + mKappaPlot;   
        mAvgMu = mAvgMu + mMeanPlot;
    end
    
    for t = 1:TMAX
        mAvgKappa(:,t) = mAvgKappa(:,t)./sum(mAvgKappa(:,t));   
        mAvgMu(:,t)    = mAvgMu(:,t)./sum(mAvgMu(:,t));
    end
  
    filename_kappa = strcat('Kappa_Decible_',num2str(numDecible),'_config_',num2str(numConfiguration),'_subject_',num2str(subj));
    savepath_1 = strcat(figure_directory_2,filename_kappa); 
    saveas(h,savepath_1);
    close all;


end

