clear all; close all; fclose all; clc

withEye                 = abs(menu('With Eye Data?', {'Yes','No'})-2);
prepost                 = menu('Conditioning?', {'Pre','Post','Conditioning'});
whichone                = menu('Which measure to use?', {'dPrime','Accuracy','RT','Dist','Efficiency','Key'});
prepostVec              = {'Precond','Postcond'};

zcrit1 = 2.5;  %but also see: https://www.nature.com/articles/s41598-018-34252-7 %%%https://core.ac.uk/download/pdf/85215435.pdf %%%%%http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.454.3890&rep=rep1&type=pdf%%%http://d-scholarship.pitt.edu/7948/1/Seo.pdf
zcrit2 = 9*58/38*0.1;  %%%%mistake in computation of PPD

prev_highreward = 1;

%% Load DATA

    rootdir= 'W:\apoores\Data\fromRoman\Jessica\Task_relevant';
    
    start_path = (rootdir);
    myPath = (start_path);
    folders = dir(myPath);
    experiment_name = 'Exp10';
    
    % JESSICA'S SUBJECTS
    subjects ={'51339677','74892169','174636117','308371706','533440903','651170259','770883157','2308149','29628872','244855064','89397817','997506577','935749766','274642591','68220649',...
        '460365483','361836210','723249755','759768504','431826010', '198895673','326983557',...
        '836370445','894073080','404159757','30063240','466420224','304361840','893517963','996801295','146624172','883218578','553187045','848451010','824648861','594360327'}; %%%% should be out headphone wrong '635418447',,'691885551' poor performance

%% Predefine Variables
s.subAccu = nan(length(subjects),1);
s.subNeut = nan(length(subjects),1);
s.subBCvH = nan(length(subjects),1);
s.subBCvL = nan(length(subjects),1);
s.subSPvH = nan(length(subjects),1);
s.subSPvL = nan(length(subjects),1);

s.Nremovedtrials= nan(length(subjects),1);
s.Nalltrials= nan(length(subjects),1);
s.condNTrl =nan(length(subjects),2,5);
s.vecAccu = cell(length(subjects),5);
s.vecAccu0 = cell(length(subjects),5);
s.vecRT0 = cell(length(subjects),5);
s.vecCondRew0 = cell(length(subjects),5);
s.vecCondMod0  = cell(length(subjects),5);
s.vecRT = cell(length(subjects),5);

for i=1:length(subjects)
   
    if withEye, load([myPath filesep subjects{i} filesep subjects{i} '__allDataEyE.mat']); vecDist = [];
    else        load([myPath filesep subjects{i} filesep subjects{i} '__allData.mat']);
    end
    vecAccu =[];  vecSoPi =[]; vecGabOALL = [];
    vecCond =[];  vecBoxC =[]; vecBoxS = [];
    vecSoVo =[];  vecValu =[];
    vecGabO =[];  vecRelP =[]; vecRT=[];
    vccboxXLoc= []; vccgaborXLoc= [];
    vectonset= [];
    vecError =[];
    keyPressed = [];

    rel_blocks =2:length(block);
    
    for j=rel_blocks % gather information into vectors
        if ~isempty(block(j).trials)
            
            RelTrl = [block(j).trials.BlPrePost] == prepost;% & [block(j).trials.RT] > 0.1 & [block(j).trials.RT] < 2.25;%&[block(j).trials.BlNumber] ~= 1;
            tempSt = sum([[block(j).trials.SL]; [block(j).trials.BL]]);
            
            vecAccu = [vecAccu block(j).trials([block(j).trials.error]==0&RelTrl).Accuracy]; % Accuracy keyPressed
            vecRT_zcrit = []; % Jessica: applied same logic as in pupil data
            vecRT_zcrit   = [vecRT_zcrit   block(j).trials([block(j).trials.error]==0& ([block(j).trials.BlPrePost] == 1 | [block(j).trials.BlPrePost] == 2)).RT]; %%%%Gabor trials RT to remove long RTs
            vecRT   = [vecRT   block(j).trials([block(j).trials.error]==0&RelTrl).RT];
            vecSoPi = [vecSoPi block(j).trials([block(j).trials.error]==0&RelTrl).SP];
            vecBoxC = [vecBoxC block(j).trials([block(j).trials.error]==0&RelTrl).BC];
            vecCond = [vecCond block(j).trials([block(j).trials.error]==0&RelTrl).condition];
            vecValu = [vecValu block(j).trials([block(j).trials.error]==0&RelTrl).StimValue];
            vecBoxS = vecGabO;
            if withEye, vecDist = [vecDist block(j).trials([block(j).trials.error]==0&RelTrl).dist];end
            vecGabOALL = [vecGabOALL block(j).trials(RelTrl).gaborOri];
            vecError = [vecError block(j).trials(RelTrl).error];
            vecGabO = [vecGabO block(j).trials([block(j).trials.error]==0&RelTrl).gaborOri];
            
            fprintf('Block:%d Rel:%d\n',j,sum(RelTrl));
        end
    end
    
    if withEye
        trls = find((vecRT>= mean(vecRT_zcrit)- zcrit1*std(vecRT_zcrit) & vecRT<mean(vecRT_zcrit)+zcrit1*std(vecRT_zcrit)) & vecDist<zcrit2);
    else
        trls = find((vecRT>= mean(vecRT_zcrit)- zcrit1*std(vecRT_zcrit) & vecRT<mean(vecRT_zcrit)+zcrit1*std(vecRT_zcrit)));
    end
    s.Nremovedtrials(i)=length(vecRT)-length(trls); % how many trials are removed
    s.Nalltrials(i)=length(vecRT); % how many trials are removed
   
    if withEye,
        fprintf('%d | %s : %d %d\n',length(subjects)-i,subjects{i},...
        length(find(abs(vecDist)>zcrit2)),length(find(abs(vecDist)<zcrit2)));
    end
    
    vecAccu = vecAccu(trls);
    vecRT   = vecRT(trls);
    if withEye,vecDist = vecDist(trls);end
    vecSoPi = vecSoPi(trls);
    vecBoxC = vecBoxC(trls);
    vecCond = vecCond(trls);
    vecValu = vecValu(trls);
    vecGabO = vecGabO(trls);
    
    % to remove outliers
    outlier = abs(zscore(vecRT))>zcrit1; % to remove slow RTs use>4;
    
    % fill with nanmean values of desired conditions across subjects
    s.subAccu(i)  = nanmean(vecAccu);                          % accuracy
    s.subRT(i)    = nanmean(vecRT);                            % RT
    % s.subminRT(i)    = nanmin(vecRT); % RT
      
    hits       = nanmean(vecAccu(vecGabO== 1));
    fA         = nanmean( 1 - vecAccu(vecGabO==-1));
    N = sum(vecGabO==-1); % for correction of h and fA when they are 1 or 0 see http://www.kangleelab.com/sdt-d-prime-calculation---other-tips.html
    s.subDPri(i)  = Dprime2(hits,fA,N);           % dPrime
    s.subRewV(i,:)= inf.rewardValueB;             % Value of cues (mode)
    s.subRewS(i,:)= inf.rewardValueS;
    s.subleanV{i}= [inf.learn.Visual];
    s.age{i}= [inf.age];
    if isfield(inf.learn,'Sound')
        s.subleanS{i}= [inf.learn.Sound];
    else
        s.subleanS{i}= 1;
    end
    s.subPSE(i)= unique([block(2).trials.PSE]);
    %% Process in dPrime
    selection{1} = (vecCond==0);                            % SubNeut
    
    selection{2} = (vecCond==1 & vecBoxC==1);               % SubBC_O
    selection{3} = (vecCond==1 & vecBoxC==2);               % SubBC_M
    selection{4} = (vecCond==1 & vecValu==1);               % SubBCvL
    selection{5} = (vecCond==1 & vecValu==2);               % SubBCvH
    
    selection{6} = (vecCond==2 & vecSoPi==1);               % SubSP_L
    selection{7} = (vecCond==2 & vecSoPi==2);               % SubSP_H
    selection{8} = (vecCond==2 & vecValu==1);               % SubSPvL
    selection{9} = (vecCond==2 & vecValu==2);               % SubSPvH
    
    selection{10} = (vecCond==1 );                          % SubV
    selection{11} = (vecCond==2 );                          % SubA
    
    forProcess = {'subNeut';...
        'subBC_O';'subBC_M';'subBCvL';'subBCvH';...
        'subSP_L';'subSP_H';'subSPvL';'subSPvH';...
        'subV';'subA'};
    
    col =1;
    for step = 1:length(forProcess)
        
        mytrls = find(selection{step});
        ntrl1 = length(mytrls);
        
        if prev_highreward
            myrealtrals =[];
            if ~isempty(mytrls)
                tcount=0;
                for t =1:length(mytrls)
                    if mytrls(t)~=1 && vecValu(mytrls(t)-1)==1 && vecCond(mytrls(t)-1)== 1 %&& vecCond(mytrls(t))== 1 %&& vecValu(mytrls(t)-2)==1 && vecCond(mytrls(t)-2)==1%vecAccu(mytrls(t)-1)==1
%                         if vecValu(mytrls(t)-2)==1 && vecCond(mytrls(t)-2)== 1
                        tcount=tcount+1;
                        myrealtrals(tcount) = mytrls(t);
%                         end
                    end
                end
            end
            mytrls=  myrealtrals;
        end
        
        hits             =   nanmean(vecAccu(mytrls(vecGabO(mytrls)== 1)));
        fA               =   nanmean(1-nanmean(vecAccu(mytrls(vecGabO(mytrls)== -1))));
       
        N = 50;
        if whichone ==1
            s.(forProcess{step})(i) = Dprime2(hits,fA,N); % dPrime
        elseif whichone ==2
            s.(forProcess{step})(i) = nanmean(vecAccu(mytrls));
        elseif whichone ==3
            s.(forProcess{step})(i) = nanmean(vecRT(mytrls));
        elseif whichone ==4
            s.(forProcess{step})(i) = nanmean(vecDist(mytrls));
        elseif whichone ==5
            s.(forProcess{step})(i) = nanmean(vecAccu(mytrls))./nanmean(vecRT(mytrls)); % efficiency score, see https://www.journalofcognition.org/articles/10.5334/joc.6/;
        elseif whichone ==6
            s.(forProcess{step})(i) = nanmean(vecAccu(mytrls)); 
        end
    end

end

%%
close all hidden
Pp = [1:length(subjects)];
for i=1:length(subjects)
    sl= [s.subleanS{i}];
    vl= [s.subleanV{i}];
    vlearn(i) = mean([vl(1:2:end)]);
    slearn(i) = mean([sl(1:2:end)]);
end
Pp(s.subAccu>0.95 | s.subAccu<0.55 | (slearn+vlearn)'/2<0.55 | (s.Nremovedtrials)./s.Nalltrials>0.2)=[];
Pp(5)=[];
rel_blocks = [2];
if betweenCond == 1
    %Do nothing
elseif betweenCond == 2
    Pp = Pp(s.subRewS(Pp,2)==2);
elseif betweenCond == 3
    Pp = Pp(s.subRewS(Pp,2)==1);
elseif betweenCond == 4
    Pp = Pp(s.subRewV(Pp,2)==2);
elseif betweenCond == 5
    Pp = Pp(s.subRewV(Pp,2)==1);
end

%% Figure 100 change of effect size across trials for visual
%
% figure(777), axis square, hold on;
% effect_mean = (nanmean(s.effect_ColorValue(Pp,:),1));
% eff_sem=(nanstd(s.effect_ColorValue(Pp,:),0,1)./sqrt(length(Pp)));
% eff_sem=1.*eff_sem; %%%% for confiedence interval (now is 90% for 95% should multiply by 1.96)
% 
% y1=effect_mean-eff_sem; y2=effect_mean+eff_sem;
% y2(isnan(y1))=[];y1(isnan(y1))=[];
% c1 = plot(1:size(effect_mean,2), nanmean(effect_mean,1),'g','linewidth',3)
% h=patch([1:length(y1) fliplr(1:length(y1))], [y2  fliplr(y1)],[0.5 0.8 0.8],'FaceAlpha',.5,'EdgeAlpha',.3)
% h.DisplayName = 'Visual';
% plot([1:size(effect_mean,2)],zeros(1,size(effect_mean,2)),':k','linewidth',2)
% axis([-1 80 -0.15 0.15])
% xlabel('Trials')
% ylabel('Performance Difference High reward minus Low Reward')
% for t=1:size(s.effect_SoundValue(Pp,:),2)
%     [h pp(t)]=ttest(s.effect_SoundValue(Pp,t),0);
%     if pp(t)<0.05
%         plot(t,-0.5,'sk','MarkerFaceColor',[0 0 0])
%     end
% end
% 
% effect_mean = (nanmean(s.effect_SoundValue(Pp,:),1));
% eff_sem=(nanstd(s.effect_SoundValue(Pp,:),0,1)./sqrt(length(Pp)));
% eff_sem=1.*eff_sem; %%%% for confiedence interval (now is 90% for 95% should multiply by 1.96)
% 
% y1=effect_mean-eff_sem; y2=effect_mean+eff_sem;
% y2(isnan(y1))=[];y1(isnan(y1))=[];
% c2 = plot(1:size(effect_mean,2), nanmean(effect_mean,1),'m','linewidth',3)
% h=patch([1:length(y1) fliplr(1:length(y1))], [y2  fliplr(y1)],[0.8 0.5 0.8],'FaceAlpha',.5,'EdgeAlpha',.3);
% h.DisplayName = 'Sound';
% plot([1:size(effect_mean,2)],zeros(1,size(effect_mean,2)),':k','linewidth',2)
% axis([-1 80 -0.15 0.15])
% 
% xlabel('ntrials')
% ylabel('Performance Difference High reward minus Low Reward')
% for t=1:size(s.effect_SoundValue(Pp,:),2)
%     [h pp(t)]=ttest(s.effect_SoundValue(Pp,t),0);
%     if pp(t)<0.05
%         plot(t,-0.5,'sk','MarkerFaceColor',[1 0 0])
%     end
% end

% %% Figure 1 Sound Pitch High-Low value
% figure(1);
% bar(1:3,[nanmean(s.subNeut(Pp)) nanmean(s.subSP_H(Pp))  nanmean(s.subSP_L(Pp))],'r'), hold on
% errorbar(1:3,[nanmean(s.subNeut(Pp)) nanmean(s.subSP_H(Pp))  nanmean(s.subSP_L(Pp))],...
%     [nanstd(s.subNeut(Pp)) nanstd(s.subSP_H(Pp)) nanstd(s.subSP_L(Pp))]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none');
% ax = gca;
% ax.XTickLabel = {'neutral','HighPitch','LowPitch'};
% if whichone==1
%     ax.YLim = scalesel1;
%     ylabel('DPrime' , 'FontSize', 14);
% elseif whichone==2
%     ax.YLim = scalesel2;
%     ylabel('Accuracy' , 'FontSize', 14);
% elseif whichone==3
%     ax.YLim = scalesel3;
%     ylabel('RT' , 'FontSize', 14);
% elseif whichone==5
%     ax.YLim = scalesel5;
%     ylabel('Efficiency' , 'FontSize', 14);
% elseif whichone==6
%     ax.YLim = scalesel6;
%     ylabel('Key (1-UP, 2-Down)' , 'FontSize', 14);
% end
% [h p]= ttest(s.subSP_H(Pp),s.subSP_L(Pp));
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', P= ' num2str((p))])


%% Figure 3 Box Color Orange-Magenta value
% figure(3);
% bar(1:3,[nanmean(s.subNeut(Pp)) nanmean(s.subBC_M(Pp))  nanmean(s.subBC_O(Pp))],'r'), hold on
% errorbar(1:3,[nanmean(s.subNeut(Pp)) nanmean(s.subBC_M(Pp))  nanmean(s.subBC_O(Pp))],...
%     [nanstd(s.subNeut(Pp)) nanstd(s.subBC_M(Pp)) nanstd(s.subBC_O(Pp))]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none')
% ax = gca;
% ax.XTickLabel = {'neutral','Magenta','Orange'};
% if whichone==1
%     ax.YLim = scalesel1;
%     ylabel('DPrime' , 'FontSize', 14);
% elseif whichone==2
%     ax.YLim = scalesel2;
%     ylabel('Accuracy' , 'FontSize', 14);
% elseif whichone==3
%     ax.YLim = scalesel3;
%     ylabel('RT' , 'FontSize', 14);
% elseif whichone==5
%     ax.YLim = scalesel5;
%     ylabel('Efficiency' , 'FontSize', 14);
% elseif whichone==6
%     ax.YLim = scalesel6;
%     ylabel('Key (1-UP, 2-Down)' , 'FontSize', 14);
% end
% [h p]= ttest( s.subBC_M(Pp),s.subBC_O(Pp));
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', P= ' num2str(round(p,2))])

%% Figure 2 Sound
figure(17);
bar([1 3 4 6 7],[nanmean(s.subNeut(Pp)) nanmean(s.subBCvH(Pp)) nanmean(s.subBCvL(Pp)) nanmean(s.subSPvH(Pp))...
    nanmean(s.subSPvL(Pp))],'w'), axis square, hold on
errorbar([1 3 4 6 7],[nanmean(s.subNeut(Pp)) nanmean(s.subBCvH(Pp)) nanmean(s.subBCvL(Pp)) nanmean(s.subSPvH(Pp))...
    nanmean(s.subSPvL(Pp))],...
    [nanstd(s.subNeut(Pp)) nanstd(s.subBCvH(Pp)) nanstd(s.subBCvL(Pp)) nanstd(s.subSPvH(Pp))...
    nanstd(s.subSPvL(Pp))]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none')
ax = gca;
ax.XTickLabel = {'Neut','VisHighVal','VisLowVal','AudioHighVal','AudioLowVal'};
ax = gca;
ax.YLim = [.7 0.85];
title('Effect of Reward value and modality', 'FontSize', 16);
[h p1]= ttest(s.subBCvH(Pp),s.subBCvL(Pp));
[h p2]= ttest(s.subSPvH(Pp),s.subSPvL(Pp));
title([ experiment_name ', Blocks:' prepostVec{prepost} ' '...
    ', N= ' num2str(length(Pp)) ', visual P= ' num2str(round(p1,3)) ', auditory P= ' num2str(round(p2,3))]);
% print ('-f2' ,'-dpsc2',[rootdir_fig 'valsounds.ps'])
%%
%% Figure 2 Sound
% figure(18);
% bar([1 2 4 5 ],[ nanmean(s.subBCvH(Pp)-s.subNeut(Pp)) nanmean(s.subBCvL(Pp)-s.subNeut(Pp)) nanmean(s.subSPvH(Pp)-s.subNeut(Pp))...
%     nanmean(s.subSPvL(Pp)-s.subNeut(Pp))],'w'), axis square, hold on
% errorbar([1 3 4 6 7],[nanmean(s.subNeut(Pp)) nanmean(s.subBCvH(Pp)) nanmean(s.subBCvL(Pp)) nanmean(s.subSPvH(Pp))...
%     nanmean(s.subSPvL(Pp))],...
%     [nanstd(s.subNeut(Pp)) nanstd(s.subBCvH(Pp)) nanstd(s.subBCvL(Pp)) nanstd(s.subSPvH(Pp))...
%     nanstd(s.subSPvL(Pp))]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none')
% ax = gca;
% ax.XTickLabel = {'Neut','VisHighVal','VisLowVal','AudioHighVal','AudioLowVal'};
% if whichone==1
%     ax.YLim = scalesel1;
%     ylabel('DPrime' , 'FontSize', 14);
% elseif whichone==2
%     ax.YLim = scalesel2;
%     ylabel('Accuracy' , 'FontSize', 14);
% elseif whichone==3
%     ax.YLim = scalesel3;
%     ylabel('RT' , 'FontSize', 14);
% elseif whichone==5
%     ax.YLim = scalesel5;
%     ylabel('Efficiency' , 'FontSize', 14);
% elseif whichone==6
%     ax.YLim = scalesel6;
%     ylabel('Key (1-UP, 2-Down)' , 'FontSize', 14);
% end
% ax = gca;
% % ax.YLim = [1.4 1.6];
% title('Effect of Reward value and modality', 'FontSize', 16);
% [h p1]= ttest(s.subBCvH(Pp),s.subBCvL(Pp));
% [h p2]= ttest(s.subSPvH(Pp),s.subSPvL(Pp));
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', visual P= ' num2str(round(p1,3)) ', auditory P= ' num2str(round(p2,3))]);
% % print ('-f2' ,'-dpsc2',[rootdir_fig 'valsounds.ps'])
% %% Figure 2 Sound
% figure(170);
% boxplot([(s.subAccu(Pp)) (s.subBCvH(Pp)) (s.subBCvL(Pp)) (s.subSPvH(Pp))...
%     (s.subSPvL(Pp))]), hold on
% % plot(1:5,[(s.subNeut(Pp)) (s.subBCvH(Pp)) (s.subBCvL(Pp)) (s.subSPvH(Pp))...
% %     (s.subSPvL(Pp))],'ok'), hold on
%
% % ax = gca;
% % ax.XTickLabel = {'Neut','VisHighVal','VisLowVal','AudioHighVal','AudioLowVal'};
% % if whichone==1
% %     ax.YLim = scalesel1;
% %     ylabel('DPrime' , 'FontSize', 14);
% % elseif whichone==2
% %     ax.YLim = scalesel2;
% %     ylabel('Accuracy' , 'FontSize', 14);
% % elseif whichone==3
% %     ax.YLim = scalesel3;
% %     ylabel('RT' , 'FontSize', 14);
% % elseif whichone==5
% %     ax.YLim = scalesel5;
% %     ylabel('Efficiency' , 'FontSize', 14);
% % elseif whichone==6
% %     ax.YLim = scalesel6;
% %     ylabel('Key (1-UP, 2-Down)' , 'FontSize', 14);
% % end
% % ax = gca;
% % ax.YLim = [1.4 1.6];
% %ylim([0.5 1])
% title('Effect of Reward value and modality', 'FontSize', 16);
% [h p1]= ttest(s.subBCvH(Pp),s.subBCvL(Pp));
% [h p2]= ttest(s.subSPvH(Pp),s.subSPvL(Pp));
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', visual P= ' num2str(round(p1,2)) ', auditory P= ' num2str(round(p2,2))]);

% %% Figure 2 Sound
% figure(180);
% boxplot([ (s.subBCvH(Pp))-(s.subNeut(Pp)) (s.subBCvL(Pp))-(s.subNeut(Pp)) (s.subSPvH(Pp))-(s.subNeut(Pp))...
%     (s.subSPvL(Pp))-(s.subNeut(Pp))]), hold on

%%
HV= s.subBCvH(Pp)';
LV= s.subBCvL(Pp)';
HA= s.subSPvH(Pp)';
LA= s.subSPvL(Pp)';
N= s.subNeut(Pp)';
anovdata_rew_modality_names = {'HV';'LV';'HA';'LA'};
anovdata_rew_modality = [HV;LV;HA;LA]';

bimodncond=length(anovdata_rew_modality_names);
subfactor=cell(bimodncond,1);

for  i=1:bimodncond
    subfactor{i}= ['F' num2str(i)];
end

t = array2table(anovdata_rew_modality,'VariableNames',subfactor);
%factorNames = {'modality','reward','latencybin','part'};
factorNames = {'modality','reward'}; % for each part


within = table({'V';'V';'A';'A'},... %modality: visual, auditory, bimodal
    {'H';'L';'H';'L'},... %reward: high or low
    'VariableNames',factorNames);
%
% fit the repeated measures model
rm = fitrm(t,['F1-F' num2str(bimodncond) '~1'],'WithinDesign',within);

% run my repeated measures anova here
%[ranovatbl] = ranova(rm, 'WithinModel','modality*reward*latencybin*part')
[ranovatbl] = ranova(rm, 'WithinModel','modality*reward')
%%
%% To include all trials, solution based on anovan, not sure how to implement in fitrm
% HV= s.subBCvH(Pp)';
% LV= s.subBCvL(Pp)';
% HA= s.subSPvH(Pp)';
% LA= s.subSPvL(Pp)';
% N= s.subNeut(Pp)';
% sel = 1:length(Pp);
% nsel=length(sel);
% anovdata = [HV;LV;HA;LA]';
% ncond=size(anovdata,2);
%
% rewconds={'H';'L';'H';'L'}';
% modalityconds={'Visual';'Visual';'Auditory';'Auditory'}';
%
% rewfactor = repmat(rewconds,nsel,1); %%%%
% modalityfactor = repmat(modalityconds,nsel,1); %%%%
%
% subj = repmat(sel,[1,ncond]);
%
% %%%%%now reshape them: column vector
% anovdata = anovdata(:);
% rewfactor = rewfactor(:); %%%%
% modalityfactor = modalityfactor(:); %%%%
% subj = subj(:); %%%%;
%
% [h p]= anovan(anovdata,{subj rewfactor  modalityfactor  },'model',[1 0 0;0 1 0;0 0 1;0 1 1],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  }); %%%%% repeated measure ANOVA
%%
% pre=load('pre_beh5_jessica.mat','s','Pp');
% HV_pre= pre.s.subBCvH(Pp)';
% LV_pre= pre.s.subBCvL(Pp)';
% HA_pre= pre.s.subSPvH(Pp)';
% LA_pre= pre.s.subSPvL(Pp)';
% N_pre= pre.s.subNeut(Pp)';
%
% % HV_pre= zeros(length(Pp),1);
% % LV_pre= zeros(length(Pp),1);
% % HA_pre= zeros(length(Pp),1);
% % LA_pre= zeros(length(Pp),1);
% % N_pre= zeros(length(Pp),1);
%
% maxn = max(max(cellfun(@length,s.vecAccu)));
% HV = nan(length(Pp),maxn);
% LV = nan(length(Pp),maxn);
% HA = nan(length(Pp),maxn);
% LA = nan(length(Pp),maxn);
% N = nan(length(Pp),maxn);
% subj = nan(maxn,length(Pp),4);
% rewconds={'H';'L';'H';'L'}';
% modalityconds={'Visual';'Visual';'Auditory';'Auditory'}';
% binconds={'1st';'2ns'}';
% prevtrial = {'high';'low';'neutr'};
% rewfactor = cell(maxn,length(Pp),length(rewconds));
% modalityfactor = cell(maxn,length(Pp),length(rewconds));
% binfactor = cell(maxn,length(Pp),length(binconds));
% prevtrialvisual = cell(maxn,length(Pp),length(binconds));
% prevtrialauditory = cell(maxn,length(Pp),length(binconds));
%
% for p=1:length(Pp)
%     HV(p,1:length(s.vecAccu{Pp(p),3}))= s.vecAccu{Pp(p),3}-HV_pre(p);
%     LV(p,1:length(s.vecAccu{Pp(p),2}))= s.vecAccu{Pp(p),2}-LV_pre(p);
%     HA(p,1:length(s.vecAccu{Pp(p),5}))= s.vecAccu{Pp(p),5}-HA_pre(p);
%     LA(p,1:length(s.vecAccu{Pp(p),4}))=s.vecAccu{Pp(p),4}-LA_pre(p);
%     N(p,1:length(s.vecAccu{Pp(p),1}))= s.vecAccu{Pp(p),1}-N_pre(p);
%     subj(1:maxn,p,1:4)= p;
%     for j =1:length(rewconds)
%         rewfactor(1:maxn,p,j) = [ repmat(rewconds(j),maxn,1)]; %%%%
%         modalityfactor(1:maxn,p,j) = [ repmat(modalityconds(j),maxn,1)]; %%%%
%         binfactor(1:maxn,p,j) = [ repmat(binconds(1),floor(maxn/2),1) ; repmat(binconds(2),ceil(maxn/2),1)]; %%%%
%         for n =1:maxn
%             if n== 1
%                 prevtrialvisual(n,p,j) = prevtrial(3); %%%%
%             elseif n<= length(s.vecAccu{Pp(p),3})
%                 if s.vecCondRew0{Pp(p),3}(n)  == 2
%                     prevtrialvisual (n,p,j) = prevtrial(1); %%%%
%                 elseif s.vecCondRew0{Pp(p),3}(n)  == 1
%                     prevtrialvisual (n,p,j) = prevtrial(2); %%%%
%                 else%if s.vecCondRew0{Pp(p),3}(n)  == 0
%                     prevtrialvisual (n,p,j) = prevtrial(3); %%%%
%                 end
%             elseif n> length(s.vecCondRew0{Pp(p),3})
%                 prevtrialvisual (n,p,j) = prevtrial(3); %%%%
%             end
%              if n== 1
%                 prevtrialauditory(n,p,j) = prevtrial(3); %%%%
%              elseif n<= length(s.vecCondRew0{Pp(p),5})
%                 if s.vecCondRew0{Pp(p),5}(n) == 2
%                     prevtrialauditory(n,p,j) =  prevtrial(1); %%%%
%                 elseif s.vecCondRew0{Pp(p),5}(n)  == 1
%                     prevtrialauditory(n,p,j) = prevtrial(2); %%%%
%                 else%if s.vecCondRew0{Pp(p),5}(n)  == 0
%                     prevtrialauditory(n,p,j) = prevtrial(3); %%%%
%                 end
%               elseif n> length(s.vecCondRew0{Pp(p),5})
%                 prevtrialauditory (n,p,j) = prevtrial(3); %%%%
%             end
%         end
%
%     end
% end
% anovdata = cat(3,HV',LV', HA', LA');
%
% %%%%%now reshape them: column vector
%
% anovdata = anovdata(:);
% subj = subj(:); %%%%;
% rewfactor = rewfactor(:); %%%%
% modalityfactor = modalityfactor(:); %%%%
% binfactor = binfactor(:); %%%%
% prevtrialauditory = prevtrialauditory(:); %%%%
% prevtrialvisual = prevtrialvisual(:); %%%%
% [h p]= anovan(anovdata,{subj rewfactor  modalityfactor binfactor  },'model',[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;0 1 1 0;0 1 0 1;0 0 1 1;0 1 1 1],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'binfactor' }); %%%%% repeated measure ANOVA
%
% %%%%%%%%%%%%%this one should be worked out
% [h p]= anovan(anovdata,{subj rewfactor  modalityfactor prevtrialvisual prevtrialauditory },'model',[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;0 1 1 0 0;0 1 0 1 0;...
%     0 1 0 0 1;0 0 1 1 0;0 0 1 0 1;0 1 1 1 1;0 1 1 0 1;0 1 1 1 0],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'previousvisual'  'previousauditory' }); %%%%% repeated measure ANOVA
%
% %%%%%%%%%%%%%this one should be worked out
% [h p]= anovan(anovdata,{subj rewfactor  modalityfactor prevtrialvisual prevtrialauditory },'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'previousvisual'  'previousauditory' }); %%%%% repeated measure ANOVA
%
% [h p]= anovan(anovdata,{subj rewfactor  modalityfactor  },'model',[1 0 0;0 1 0;0 0 1;0 1 1],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  }); %%%%% repeated measure ANOVA
%%

%%
% figure(181);
% boxplot([HV-LV;HA-LA]')
% [h p1]= ttest(HV,LV);
% [h p2]= ttest(HA,LA);
% title([ experiment_name  ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', visual P= ' num2str(round(p1,2)) ', auditory P= ' num2str(round(p2,2))]);

%%
pre=load('pre_beh5_jessica.mat','s','Pp');
HV= s.subBCvH(Pp)';
LV= s.subBCvL(Pp)';
HA= s.subSPvH(Pp)';
LA= s.subSPvL(Pp)';
N= s.subNeut(Pp)';
HV_pre= pre.s.subBCvH(Pp)';
LV_pre= pre.s.subBCvL(Pp)';
HA_pre= pre.s.subSPvH(Pp)';
LA_pre= pre.s.subSPvL(Pp)';
N_pre= pre.s.subNeut(Pp)';
anovdata_rew_modality_names = {'HV';'LV';'HA';'LA';'HV_pre';'LV_pre';'HA_pre';'LA_pre'};
anovdata_rew_modality = [HV;LV;HA;LA;HV_pre;LV_pre;HA_pre;LA_pre]';

bimodncond=length(anovdata_rew_modality_names);
subfactor=cell(bimodncond,1);

for  i=1:bimodncond
    subfactor{i}= ['F' num2str(i)];
end

t = array2table([anovdata_rew_modality (1:length(Pp))'],'VariableNames',[subfactor;{'Sujs'}]);
%factorNames = {'modality','reward','latencybin','part'};
factorNames = {'modality','reward','pre_post'}; % for each part


within = table({'V';'V';'A';'A';'V';'V';'A';'A'},... %modality: visual, auditory, bimodal
    {'H';'L';'H';'L';'H';'L';'H';'L'},... %reward: high or low
    {'Post';'Post';'Post';'Post';'Pre';'Pre';'Pre';'Pre'},... %phase: pre or post
    'VariableNames',factorNames);
%
% fit the repeated measures model
rm = fitrm(t,['F1-F' num2str(bimodncond) '~1'],'WithinDesign',within);

% run my repeated measures anova here
%[ranovatbl] = ranova(rm, 'WithinModel','modality*reward*latencybin*part')
[ranovatbl] = ranova(rm, 'WithinModel','modality*reward*pre_post')


%% linear model version based on averaged data: see line 1817 for linear model based on all data
% pre=load('pre_beh5_jessica.mat','s','Pp');
% HV= s.subBCvH(Pp)';
% LV= s.subBCvL(Pp)';
% HA= s.subSPvH(Pp)';
% LA= s.subSPvL(Pp)';
% N= s.subNeut(Pp)';
% HV_pre= pre.s.subBCvH(Pp)';
% LV_pre= pre.s.subBCvL(Pp)';
% HA_pre= pre.s.subSPvH(Pp)';
% LA_pre= pre.s.subSPvL(Pp)';
% N_pre= pre.s.subNeut(Pp)';
% anovdata_rew_modality = [HV;LV;HA;LA;HV_pre;LV_pre;HA_pre;LA_pre]';
%
% rewconds={'H';'L';'H';'L';'H';'L';'H';'L'}';
% modalityconds={'Visual';'Visual';'Auditory';'Auditory';'Visual';'Visual';'Auditory';'Auditory'}';
% pconds={'Post';'Post';'Post';'Post';'Pre';'Pre';'Pre';'Pre'}';
% nsel = length(Pp);
% ncond = length(rewconds);
% rewfactor = repmat(rewconds,nsel,1); %%%%
% modalityfactor = repmat(modalityconds,nsel,1); %%%%
% pfactor = repmat(pconds,nsel,1); %%%%
%
% subj = repmat(1:length(Pp),[1,ncond]);
%
% %%%%%now reshape them: column vector
% anovdata = anovdata_rew_modality(:);
% rewfactor = rewfactor(:); %%%%
% modalityfactor = modalityfactor(:); %%%%
% pfactor = pfactor(:);
% subj = subj(:); %%%%;
% tbl2 = table(subj,modalityfactor,rewfactor,pfactor,anovdata,'VariableNames',{'subjects'  'modalityfactor','rewardfactor','phase','accuracy' });
% lm = fitlm(tbl2,'accuracy~1+rewardfactor*modalityfactor*phase')

%%
pre=load('pre_beh5_jessica.mat','s','Pp');

HV= s.subBCvH(Pp)';
LV= s.subBCvL(Pp)';
HA= s.subSPvH(Pp)';
LA= s.subSPvL(Pp)';
N= s.subNeut(Pp)';


HV_pre= pre.s.subBCvH(Pp)';
LV_pre= pre.s.subBCvL(Pp)';
HA_pre= pre.s.subSPvH(Pp)';
LA_pre= pre.s.subSPvL(Pp)';
N_pre= pre.s.subNeut(Pp)';


post1H=load('post_beh5_jessica_afterHigh_nminus1.mat','s','Pp');
post2H=load('post_beh5_jessica_afterHigh2.mat','s','Pp');

post1L=load('post_beh5_jessica_afterLow_nminus1.mat','s','Pp');
post2L=load('post_beh5_jessica_afterLow2.mat','s','Pp');


n1h= (post1H.s.subBCvH(Pp)'+post1H.s.subBCvL(Pp)'+post1H.s.subSPvH(Pp)'+ post1H.s.subSPvL(Pp)'+post1H.s.subNeut(Pp)')/5;

n1l= (post1L.s.subBCvH(Pp)'+post1L.s.subBCvL(Pp)'+post1L.s.subSPvH(Pp)'+ post1L.s.subSPvL(Pp)'+post1L.s.subNeut(Pp)')/5;


HV1= post2H.s.subBCvH(Pp)';
LV1= post2H.s.subBCvL(Pp)';
HA1= post2H.s.subSPvH(Pp)';
LA1= post2H.s.subSPvL(Pp)';
N1= post2H.s.subNeut(Pp)';

HV2= post2L.s.subBCvH(Pp)';
LV2= post2L.s.subBCvL(Pp)';
HA2= post2L.s.subSPvH(Pp)';
LA2= post2L.s.subSPvL(Pp)';
N2= post2L.s.subNeut(Pp)';


% figure(1732), axis square, hold on;
% boxplot([((HV+HA)./2-(HV_pre+HA_pre)./2)'  ((LV+LA)./2-(LV_pre+LA_pre)./2)' (N-N_pre)'])


figure(1733), axis square, hold on;
% bar([1 ],[ nanmean((HV+HA)./2-(HV_pre+HA_pre)./2)  ],'r'), hold on
%
% bar([3],[  nanmean((LV+LA)./2-(LV_pre+LA_pre)./2) ],'b')
%
% bar([5 ],[  nanmean(N-N_pre)],'k')
errorbar([1 ],[ nanmean((HV+HA)./2-(HV_pre+HA_pre)./2)  ],...
    [nanstd((HV+HA)./2-(HV_pre+HA_pre)./2)  ]./sqrt(length(Pp)),'rd','linewidth',3), hold on
% errorbar([1.1 ],[ nanmean((HA+HA)./2-(HA_pre+HA_pre)./2)  ],...
%     [nanstd((HA+HA)./2-(HA_pre+HA_pre)./2)  ]./sqrt(length(Pp)),'mx','linewidth',3), hold on
% errorbar([1.2 ],[ nanmean((HV+HV)./2-(HV_pre+HV_pre)./2)  ],...
%     [nanstd((HV+HV)./2-(HV_pre+HV_pre)./2)  ]./sqrt(length(Pp)),'cx','linewidth',3), hold on
% errorbar([1.3 ],[ nanmean((HA1+HA1)./2-(HA_pre+HA_pre)./2)  ],...
%     [nanstd((HA1+HA1)./2-(HA_pre+HA_pre)./2)  ]./sqrt(length(Pp)),'mx','linewidth',1), hold on
% errorbar([1.4 ],[ nanmean((HA2+HA2)./2-(HA_pre+HA_pre)./2)  ],...
%     [nanstd((HA2+HA2)./2-(HA_pre+HA_pre)./2)  ]./sqrt(length(Pp)),'mx','linewidth',1), hold on
% errorbar([1.2 ],[ nanmean((HV+HV)./2-(HV_pre+HV_pre)./2)  ],...
%     [nanstd((HV+HV)./2-(HV_pre+HV_pre)./2)  ]./sqrt(length(Pp)),'cx','linewidth',3), hold on
% errorbar([1.4 ],[ nanmean((HV1+HV1)./2-(HV_pre+HV_pre)./2)  ],...
%     [nanstd((HV1+HV1)./2-(HV_pre+HV_pre)./2)  ]./sqrt(length(Pp)),'cx','linewidth',1), hold on
% errorbar([1.5 ],[ nanmean((HV2+HV2)./2-(HV_pre+HV_pre)./2)  ],...
%     [nanstd((HV2+HV2)./2-(HV_pre+HV_pre)./2)  ]./sqrt(length(Pp)),'cx','linewidth',1), hold on
% 


errorbar([3],[  nanmean((LV+LA)./2-(LV_pre+LA_pre)./2) ],...
    [ nanstd((LV+LA)./2-(LV_pre+LA_pre)./2)  ]./sqrt(length(Pp)),'bd','linewidth',3)
% errorbar([3.1 ],[ nanmean(LA-LA_pre) ],...
%     [nanstd(LA-LA_pre)]./sqrt(length(Pp)),'mx','linewidth',3), hold on
% errorbar([3.4 ],[ nanmean((LA1+LA1)./2-(LA_pre+LA_pre)./2)  ],...
%     [nanstd((LA1+LA1)./2-(LA_pre+LA_pre)./2)  ]./sqrt(length(Pp)),'mx','linewidth',1), hold on
% errorbar([3.3 ],[ nanmean((LA2+LA2)./2-(LA_pre+LA_pre)./2)  ],...
%     [nanstd((LA2+LA2)./2-(LA_pre+LA_pre)./2)  ]./sqrt(length(Pp)),'mx','linewidth',1), hold on
% errorbar([3.2 ],[ nanmean(LV-LV_pre) ],...
%     [nanstd(LV-LV_pre)  ]./sqrt(length(Pp)),'cx','linewidth',3), hold on
% errorbar([3.4 ],[ nanmean((LV1+LV1)./2-(LV_pre+LV_pre)./2)  ],...
%     [nanstd((LV1+LV1)./2-(LV_pre+LV_pre)./2)  ]./sqrt(length(Pp)),'cx','linewidth',1), hold on
% errorbar([3.5 ],[ nanmean((LV2+LV2)./2-(LV_pre+LV_pre)./2)  ],...
%     [nanstd((LV2+LV2)./2-(LV_pre+LV_pre)./2)  ]./sqrt(length(Pp)),'cx','linewidth',1), hold on

errorbar([5 ],[  nanmean(N-N_pre)],...
    [ nanstd(N-N_pre)  ]./sqrt(length(Pp)),'kd','linewidth',3)

ax = gca;
ax.XTickLabel = {'High','Low','Neutral'};
ax.XTick = [ 1  3  5 ];
ax.YLim = [-.05 0.1];
plot([0 6],[0 0],':k')
[h p_hh]= ttest((HV+HA)./2-(HV_pre+HA_pre)./2)  %%%%% this is switch versus no switch
[h p_ll]= ttest((LV+LA)./2-(LV_pre+LA_pre)./2)  %%%%% this is switch versus no switch
[h p_n]= ttest((N)-(N_pre))  %%%%% this is switch versus no switch
[h p_hl]= ttest(((HV+HA)./2-(HV_pre+HA_pre)./2),((LV+LA)./2-(LV_pre+LA_pre)./2 ))
[h p_hla]= ttest(((HA+HA)./2-(HA_pre+HA_pre)./2),((LA+LA)./2-(LA_pre+LA_pre)./2 ))
%[h p_hl]= ttest((HV+HA)./2-(HV_pre+HA_pre)./2,N-N_pre)

[h p_hl]= ttest((HV+HA)./2-(LV+LA)./2,(HV_pre+HA_pre)./2-(LV_pre+LA_pre)./2 )
title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
    ', N= ' num2str(length(Pp)) ', Learning_H P= ' num2str(round(p_hh,3)) ', Learning_L P= ' num2str(round(p_ll,3))...
    ', Learning_N P= ' num2str(round(p_n,3)) ', Difference HIgh and Low P= ' num2str(round(p_hl,3))]);
%%
pre=load('pre_beh5_jessica.mat','s','Pp');

HV= s.subBCvH(Pp)';
LV= s.subBCvL(Pp)';
HA= s.subSPvH(Pp)';
LA= s.subSPvL(Pp)';
N= s.subNeut(Pp)';


HV_pre= pre.s.subBCvH(Pp)';
LV_pre= pre.s.subBCvL(Pp)';
HA_pre= pre.s.subSPvH(Pp)';
LA_pre= pre.s.subSPvL(Pp)';
N_pre= pre.s.subNeut(Pp)';


post1H=load('post_beh5_jessica_afterHigh_nminus1.mat','s','Pp');
post2H=load('post_beh5_jessica_afterHigh2.mat','s','Pp');

post1L=load('post_beh5_jessica_afterLow_nminus1.mat','s','Pp');
post2L=load('post_beh5_jessica_afterLow2.mat','s','Pp');


post1H_c=load('post_beh5_jessica_afterHigh_nminus1_correct.mat','s','Pp');
post2H_c=load('post_beh5_jessica_afterHigh2_correct.mat','s','Pp');

post1L_c=load('post_beh5_jessica_afterLow_nminus1_correct.mat','s','Pp');
post2L_c=load('post_beh5_jessica_afterLow2_correct.mat','s','Pp');

post1H_e=load('post_beh5_jessica_afterHigh_nminus1_error.mat','s','Pp');
post2H_e=load('post_beh5_jessica_afterHigh2_error.mat','s','Pp');

post1L_e=load('post_beh5_jessica_afterLow_nminus1_error.mat','s','Pp');
post2L_e=load('post_beh5_jessica_afterLow2_error.mat','s','Pp');


n1h= (post1H.s.subBCvH(Pp)'+post1H.s.subBCvL(Pp)'+post1H.s.subSPvH(Pp)'+ post1H.s.subSPvL(Pp)'+post1H.s.subNeut(Pp)')/5;

n1l= (post1L.s.subBCvH(Pp)'+post1L.s.subBCvL(Pp)'+post1L.s.subSPvH(Pp)'+ post1L.s.subSPvL(Pp)'+post1L.s.subNeut(Pp)')/5;


HV1= post2H.s.subBCvH(Pp)';
LV1= post2H.s.subBCvL(Pp)';
HA1= post2H.s.subSPvH(Pp)';
LA1= post2H.s.subSPvL(Pp)';
N1= post2H.s.subNeut(Pp)';

HV2= post2L.s.subBCvH(Pp)';
LV2= post2L.s.subBCvL(Pp)';
HA2= post2L.s.subSPvH(Pp)';
LA2= post2L.s.subSPvL(Pp)';
N2= post2L.s.subNeut(Pp)';


% figure(1732), axis square, hold on;
% boxplot([((HV+HA)./2-(HV_pre+HA_pre)./2)'  ((LV+LA)./2-(LV_pre+LA_pre)./2)' (N-N_pre)'])


figure(1743), axis square, hold on;
% bar([1 ],[ nanmean((HV+HA)./2-(HV_pre+HA_pre)./2)  ],'r'), hold on
%
% bar([3],[  nanmean((LV+LA)./2-(LV_pre+LA_pre)./2) ],'b')
%
% bar([5 ],[  nanmean(N-N_pre)],'k')
% errorbar([1 3],[ nanmean((HV+HA)./2-(HV_pre+HA_pre)./2)  nanmean((LV+LA)./2-(LV_pre+LA_pre)./2)],...
%     [nanstd((HV+HA)./2-(HV_pre+HA_pre)./2) nanstd((LV+LA)./2-(LV_pre+LA_pre)./2) ]./sqrt(length(Pp)),'r','linewidth',3), hold on
% errorbar([1.1 3.1],[ nanmean((HA+HA)./2-(HA_pre+HA_pre)./2) nanmean(LA-LA_pre)  ],...
%     [nanstd((HA+HA)./2-(HA_pre+HA_pre)./2)  nanstd(LA-LA_pre)]./sqrt(length(Pp)),'m','linewidth',3), hold on
% errorbar([1.2 3.2],[ nanmean((HV+HV)./2-(HV_pre+HV_pre)./2)  nanmean(LV-LV_pre) ],...
%     [nanstd((HV+HV)./2-(HV_pre+HV_pre)./2)  nanstd(LV-LV_pre)  ]./sqrt(length(Pp)),'c','linewidth',3), hold on

errorbar([1.3 3.3],[ nanmean((HA1+HA1)./2-(HA_pre+HA_pre)./2) nanmean((LA1+LA1)./2-(LA_pre+LA_pre)./2)  ],...
    [nanstd((HA1+HA1)./2-(HA_pre+HA_pre)./2)  nanstd((LA1+LA1)./2-(LA_pre+LA_pre)./2)  ]./sqrt(length(Pp)),'m','linewidth',1), hold on
errorbar([1.4 3.4],[ nanmean((HA2+HA2)./2-(HA_pre+HA_pre)./2)  nanmean((LA2+LA2)./2-(LA_pre+LA_pre)./2) ],...
    [nanstd((HA2+HA2)./2-(HA_pre+HA_pre)./2)  nanstd((LA2+LA2)./2-(LA_pre+LA_pre)./2)]./sqrt(length(Pp)),'m:','linewidth',1), hold on

errorbar([1.4 3.4],[ nanmean((HV1+HV1)./2-(HV_pre+HV_pre)./2) nanmean((LV1+LV1)./2-(LV_pre+LV_pre)./2)   ],...
    [nanstd((HV1+HV1)./2-(HV_pre+HV_pre)./2)  nanstd((LV1+LV1)./2-(LV_pre+LV_pre)./2) ]./sqrt(length(Pp)),'c','linewidth',1), hold on
errorbar([1.5 3.5],[ nanmean((HV2+HV2)./2-(HV_pre+HV_pre)./2) nanmean((LV2+LV2)./2-(LV_pre+LV_pre)./2) ],...
    [nanstd((HV2+HV2)./2-(HV_pre+HV_pre)./2) nanstd((LV2+LV2)./2-(LV_pre+LV_pre)./2)  ]./sqrt(length(Pp)),':c','linewidth',1), hold on

errorbar([5 ],[  nanmean(N-N_pre)],...
    [ nanstd(N-N_pre)  ]./sqrt(length(Pp)),'kd','linewidth',3)

ax = gca;
ax.XTickLabel = {'High','Low','Neutral'};
ax.XTick = [ 1  3  5 ];
ax.YLim = [-.05 0.1];
plot([0 6],[0 0],':k')
[h p_hh]= ttest((HV+HA)./2-(HV_pre+HA_pre)./2)  %%%%% this is switch versus no switch
[h p_ll]= ttest((LV+LA)./2-(LV_pre+LA_pre)./2)  %%%%% this is switch versus no switch
[h p_n]= ttest((N)-(N_pre))  %%%%% this is switch versus no switch
[h p_hl]= ttest(((HV+HA)./2-(HV_pre+HA_pre)./2),((LV+LA)./2-(LV_pre+LA_pre)./2 ))
[h p_hla]= ttest(((HA+HA)./2-(HA_pre+HA_pre)./2),((LA+LA)./2-(LA_pre+LA_pre)./2 ))
%[h p_hl]= ttest((HV+HA)./2-(HV_pre+HA_pre)./2,N-N_pre)
[h p_hla2]= ttest(((HA1+HV1)./2-(HA_pre+HV_pre)./2),((LA1+LV1)./2-(LA_pre+LV_pre)./2 ))
[h p_hla3]= ttest(((HA1)./1-(HA_pre)./1),((LA1)./1-(LA_pre)./1 ))
[h p_hla3]= ttest(((HA1+HV1)./2-(HA_pre+HV_pre)./2),((LA1+LV1)./2-(LA_pre+LV_pre)./2 ))


title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
    ', N= ' num2str(length(Pp)) ', Learning_H P= ' num2str(round(p_hh,3)) ', Learning_L P= ' num2str(round(p_ll,3))...
    ', Learning_N P= ' num2str(round(p_n,3)) ', Difference HIgh and Low P= ' num2str(round(p_hl,3))]);
%%
shoulddoRT =0;

if shoulddoRT
    pre=load('pre_beh5_jessica_RT.mat','s','Pp');
        
    post2H=load('post_beh5_jessica_afterHigh2_RT.mat','s','Pp');
    post1L=load('post_beh5_jessica_afterLow_nminus1.mat','s','Pp');
    post2L=load('post_beh5_jessica_afterLow2_RT.mat','s','Pp');
        

    post2H_c=load('post_beh5_jessica_afterHigh2_correct_RT.mat','s','Pp');
    post2L_c=load('post_beh5_jessica_afterLow2_correct_RT.mat','s','Pp');
    post2N_c=load('post_beh5_jessica_afterNeut2_correct_RT.mat','s','Pp');
    
        
    post2H_e=load('post_beh5_jessica_afterHigh2_error_RT.mat','s','Pp');
    post2L_e=load('post_beh5_jessica_afterLow2_error_RT.mat','s','Pp');
    post2N_e=load('post_beh5_jessica_afterNeut2_error_RT.mat','s','Pp');
    ylabel('Reaction Time')

    myYLim = [.7 0.9];   

else
    
    pre=load('pre_beh5_jessica.mat','s','Pp');
    post2H=load('post_beh5_jessica_afterHigh2.mat','s','Pp');
    post2L=load('post_beh5_jessica_afterLow2.mat','s','Pp');
    
    
    post2H_c=load('post_beh5_jessica_afterHigh2_correct.mat','s','Pp');
    post2L_c=load('post_beh5_jessica_afterLow2_correct.mat','s','Pp');
    post2N_c=load('post_beh5_jessica_afterNeut2_correct.mat','s','Pp');
    
    
    post2H_e=load('post_beh5_jessica_afterHigh2_error.mat','s','Pp');
    post2L_e=load('post_beh5_jessica_afterLow2_error.mat','s','Pp');
    post2N_e=load('post_beh5_jessica_afterNeut2_error.mat','s','Pp');
    ylabel('Performance Difference from Baseline')

    myYLim = [-0.1 0.1];
end


HV= s.subBCvH(Pp)';
LV= s.subBCvL(Pp)';
HA= s.subSPvH(Pp)';
LA= s.subSPvL(Pp)';
N= s.subNeut(Pp)';


HV_pre= pre.s.subBCvH(Pp)';
LV_pre= pre.s.subBCvL(Pp)';
HA_pre= pre.s.subSPvH(Pp)';
LA_pre= pre.s.subSPvL(Pp)';
N_pre= pre.s.subNeut(Pp)';

if shoulddoRT
HV_pre= zeros(length(Pp),1)';
LV_pre= zeros(length(Pp),1)';
HA_pre= zeros(length(Pp),1)';
LA_pre= zeros(length(Pp),1)';
N_pre= zeros(length(Pp),1)';
end

n1h= (post1H.s.subBCvH(Pp)'+post1H.s.subBCvL(Pp)'+post1H.s.subSPvH(Pp)'+ post1H.s.subSPvL(Pp)'+post1H.s.subNeut(Pp)')/5;

n1l= (post1L.s.subBCvH(Pp)'+post1L.s.subBCvL(Pp)'+post1L.s.subSPvH(Pp)'+ post1L.s.subSPvL(Pp)'+post1L.s.subNeut(Pp)')/5;


HV1= post2H.s.subBCvH(Pp)';
LV1= post2H.s.subBCvL(Pp)';
HA1= post2H.s.subSPvH(Pp)';
LA1= post2H.s.subSPvL(Pp)';
N1= post2H.s.subNeut(Pp)';

HV2= post2L.s.subBCvH(Pp)';
LV2= post2L.s.subBCvL(Pp)';
HA2= post2L.s.subSPvH(Pp)';
LA2= post2L.s.subSPvL(Pp)';
N2= post2L.s.subNeut(Pp)';


HV1_c= post2H_c.s.subBCvH(Pp)';
LV1_c= post2H_c.s.subBCvL(Pp)';
HA1_c= post2H_c.s.subSPvH(Pp)';
LA1_c= post2H_c.s.subSPvL(Pp)';
N1_c= post2H_c.s.subNeut(Pp)';

HV2_c= post2L_c.s.subBCvH(Pp)';
LV2_c= post2L_c.s.subBCvL(Pp)';
HA2_c= post2L_c.s.subSPvH(Pp)';
LA2_c= post2L_c.s.subSPvL(Pp)';
N2_c= post2L_c.s.subNeut(Pp)';

HVn_c= post2N_c.s.subBCvH(Pp)';
LVn_c= post2N_c.s.subBCvL(Pp)';
HAn_c= post2N_c.s.subSPvH(Pp)';
LAn_c= post2N_c.s.subSPvL(Pp)';
Nn_c= post2N_c.s.subNeut(Pp)';

HV1_e= post2H_e.s.subBCvH(Pp)';
LV1_e= post2H_e.s.subBCvL(Pp)';
HA1_e= post2H_e.s.subSPvH(Pp)';
LA1_e= post2H_e.s.subSPvL(Pp)';
N1_e= post2H_e.s.subNeut(Pp)';

HV2_e= post2L_e.s.subBCvH(Pp)';
LV2_e= post2L_e.s.subBCvL(Pp)';
HA2_e= post2L_e.s.subSPvH(Pp)';
LA2_e= post2L_e.s.subSPvL(Pp)';
N2_e= post2L_e.s.subNeut(Pp)';

HVn_e= post2N_e.s.subBCvH(Pp)';
LVn_e= post2N_e.s.subBCvL(Pp)';
HAn_e= post2N_e.s.subSPvH(Pp)';
LAn_e= post2N_e.s.subSPvL(Pp)';
Nn_e= post2N_e.s.subNeut(Pp)';

H1_c= (HV1_c+HA1_c)./2;
L1_c= (LV1_c+LA1_c)./2;
H2_c= (HV2_c+HA2_c)./2;
L2_c= (LV2_c+LA2_c)./2;
Hn_c= (HVn_c+HAn_c)./2;
Ln_c= (LVn_c+LAn_c)./2;

H1_e= (HV1_e+HA1_e)./2;
L1_e= (LV1_e+LA1_e)./2;
H2_e= (HV2_e+HA2_e)./2;
L2_e= (LV2_e+LA2_e)./2;
Hn_e= (HVn_e+HAn_e)./2;
Ln_e= (LVn_e+LAn_e)./2;


% figure(1732), axis square, hold on;
% boxplot([((HV+HA)./2-(HV_pre+HA_pre)./2)'  ((LV+LA)./2-(LV_pre+LA_pre)./2)' (N-N_pre)'])


% figure(17430), axis square, hold on;
% 
% errorbar([1.3 3.3],[ nanmean((HA1_c+HA1_c)./2-(HA_pre+HA_pre)./2) nanmean((LA1_c+LA1_c)./2-(LA_pre+LA_pre)./2)  ],...
%     [nanstd((HA1_c+HA1_c)./2-(HA_pre+HA_pre)./2)  nanstd((LA1_c+LA1_c)./2-(LA_pre+LA_pre)./2)  ]./sqrt(length(Pp)),'m','linewidth',2), hold on
% errorbar([1.5 3.4],[ nanmean((HA2_c+HA2_c)./2-(HA_pre+HA_pre)./2)  nanmean((LA2_c+LA2_c)./2-(LA_pre+LA_pre)./2) ],...
%     [nanstd((HA2_c+HA2_c)./2-(HA_pre+HA_pre)./2)  nanstd((LA2_c+LA2_c)./2-(LA_pre+LA_pre)./2)]./sqrt(length(Pp)),'m:','linewidth',2), hold on
% 
% 
% errorbar([1.5 3.4],[ nanmean((HV1_c+HV1_c)./2-(HV_pre+HV_pre)./2) nanmean((LV1_c+LV1_c)./2-(LV_pre+LV_pre)./2)   ],...
%     [nanstd((HV1_c+HV1_c)./2-(HV_pre+HV_pre)./2)  nanstd((LV1_c+LV1_c)./2-(LV_pre+LV_pre)./2) ]./sqrt(length(Pp)),'c','linewidth',2), hold on
% errorbar([1.5 3.5],[ nanmean((HV2_c+HV2_c)./2-(HV_pre+HV_pre)./2) nanmean((LV2_c+LV2_c)./2-(LV_pre+LV_pre)./2) ],...
%     [nanstd((HV2_c+HV2_c)./2-(HV_pre+HV_pre)./2) nanstd((LV2_c+LV2_c)./2-(LV_pre+LV_pre)./2)  ]./sqrt(length(Pp)),':c','linewidth',2), hold on
% 
% errorbar([30 ],[  nanmean(N1_c-N_pre)],...
%     [ nanstd(N1_c-N_pre)  ]./sqrt(length(Pp)),'kd','linewidth',2)
% 
% 
% 
% errorbar([1.3 3.3]+5,[ nanmean((HA1_e+HA1_e)./2-(HA_pre+HA_pre)./2) nanmean((LA1_e+LA1_e)./2-(LA_pre+LA_pre)./2)  ],...
%     [nanstd((HA1_e+HA1_e)./2-(HA_pre+HA_pre)./2)  nanstd((LA1_e+LA1_e)./2-(LA_pre+LA_pre)./2)  ]./sqrt(length(Pp)),'m','linewidth',3), hold on
% errorbar([1.5 3.4]+5,[ nanmean((HA2_e+HA2_e)./2-(HA_pre+HA_pre)./2)  nanmean((LA2_e+LA2_e)./2-(LA_pre+LA_pre)./2) ],...
%     [nanstd((HA2_e+HA2_e)./2-(HA_pre+HA_pre)./2)  nanstd((LA2_e+LA2_e)./2-(LA_pre+LA_pre)./2)]./sqrt(length(Pp)),'m:','linewidth',3), hold on
% 
% errorbar([1.5 3.5]+5,[ nanmean((HV1_e+HV1_e)./2-(HV_pre+HV_pre)./2) nanmean((LV1_e+LV1_e)./2-(LV_pre+LV_pre)./2)   ],...
%     [nanstd((HV1_e+HV1_e)./2-(HV_pre+HV_pre)./2)  nanstd((LV1_e+LV1_e)./2-(LV_pre+LV_pre)./2) ]./sqrt(length(Pp)),'c','linewidth',3), hold on
% errorbar([1.5 3.5]+5,[ nanmean((HV2_e+HV2_e)./2-(HV_pre+HV_pre)./2) nanmean((LV2_e+LV2_e)./2-(LV_pre+LV_pre)./2) ],...
%     [nanstd((HV2_e+HV2_e)./2-(HV_pre+HV_pre)./2) nanstd((LV2_e+LV2_e)./2-(LV_pre+LV_pre)./2)  ]./sqrt(length(Pp)),':c','linewidth',3), hold on
% 
% errorbar([30 ],[  nanmean(N1_e-N_pre)],...
%     [ nanstd(N1_e-N_pre)  ]./sqrt(length(Pp)),'kd','linewidth',3)
% 
% 
% errorbar([1.3 3.3]+10,[ nanmean((H1_c+H1_c)./2-(HA_pre+HV_pre)./2) nanmean((L1_c+L1_c)./2-(LA_pre+LV_pre)./2)  ],...
%     [nanstd((H1_c+H1_c)./2-(HA_pre+HV_pre)./2)  nanstd((L1_c+L1_c)./2-(LA_pre+LV_pre)./2)  ]./sqrt(length(Pp)),'r','linewidth',2), hold on
% errorbar([1.5 3.5]+10,[ nanmean((H2_c+H2_c)./2-(HA_pre+HV_pre)./2)  nanmean((L2_c+L2_c)./2-(LA_pre+LV_pre)./2) ],...
%     [nanstd((H2_c+H2_c)./2-(HA_pre+HV_pre)./2)  nanstd((L2_c+L2_c)./2-(LA_pre+LV_pre)./2)]./sqrt(length(Pp)),'b:','linewidth',2), hold on
% 
% errorbar([1.3 3.3]+20,[ nanmean((H1_e+H1_e)./2-(HA_pre+HV_pre)./2) nanmean((L1_e+L1_e)./2-(LA_pre+LV_pre)./2)  ],...
%     [nanstd((H1_e+H1_e)./2-(HA_pre+HV_pre)./2)  nanstd((L1_e+L1_e)./2-(LA_pre+LV_pre)./2)  ]./sqrt(length(Pp)),'k','linewidth',3), hold on
% errorbar([1.5 3.5]+20,[ nanmean((H2_e+H2_e)./2-(HA_pre+HV_pre)./2)  nanmean((L2_e+L2_e)./2-(LA_pre+LV_pre)./2) ],...
%     [nanstd((H2_e+H2_e)./2-(HA_pre+HV_pre)./2)  nanstd((L2_e+L2_e)./2-(LA_pre+LV_pre)./2)]./sqrt(length(Pp)),'b:','linewidth',3), hold on
% 
% ax = gca;
% ax.XTickLabel = {'High','Low','Neutral'};
% ax.XTick = [ 1  3  5 ];
% ax.YLim = myYLim;
% plot([0 6],[0 0],':k')
% [h p_hh]= ttest((HV+HA)./2-(HV_pre+HA_pre)./2)  %%%%% this is switch versus no switch
% [h p_ll]= ttest((LV+LA)./2-(LV_pre+LA_pre)./2)  %%%%% this is switch versus no switch
% [h p_n]= ttest((N)-(N_pre))  %%%%% this is switch versus no switch
% [h p_hl]= ttest(((HV+HA)./2-(HV_pre+HA_pre)./2),((LV+LA)./2-(LV_pre+LA_pre)./2 ))
% [h p_hla]= ttest(((HA+HA)./2-(HA_pre+HA_pre)./2),((LA+LA)./2-(LA_pre+LA_pre)./2 ))
% %[h p_hl]= ttest((HV+HA)./2-(HV_pre+HA_pre)./2,N-N_pre)
% [h p_hla2]= ttest(((HA1+HV1)./2-(HA_pre+HV_pre)./2),((LA1+LV1)./2-(LA_pre+LV_pre)./2 ))
% [h p_hla3]= ttest(((HA1)./1-(HA_pre)./1),((LA1)./1-(LA_pre)./1 ))
% [h p_hla4]= ttest(((HA1+HV1)./2-(HA_pre+HV_pre)./2)-((LA1+LV1)./2-(LA_pre+LV_pre)./2 ),((HA2+HV2)./2-(HA_pre+HV_pre)./2)-((LA2+LV2)./2-(LA_pre+LV_pre)./2 ))
% 
% 
% 
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', Learning_H P= ' num2str(round(p_hh,3)) ', Learning_L P= ' num2str(round(p_ll,3))...
%     ', Learning_N P= ' num2str(round(p_n,3)) ', Difference HIgh and Low P= ' num2str(round(p_hl,3))]);



figure(17431), hold on;

set(figure(17431), 'Position', [0 0 1024 768]);
set(gca,'FontSize',24);
errorbar([1.3 3.3],[ nanmean((H1_c+H1_c)./2-(HA_pre+HV_pre)./2) nanmean((L1_c+L1_c)./2-(LA_pre+LV_pre)./2)  ],...
    [nanstd((H1_c+H1_c)./2-(HA_pre+HV_pre)./2)  nanstd((L1_c+L1_c)./2-(LA_pre+LV_pre)./2)  ]./sqrt(length(Pp)),'color',[1 0.5 0],'marker','d','linewidth',6), hold on
errorbar([1.5 3.5],[ nanmean((H2_c+H2_c)./2-(HA_pre+HV_pre)./2)  nanmean((L2_c+L2_c)./2-(LA_pre+LV_pre)./2) ],...
    [nanstd((H2_c+H2_c)./2-(HA_pre+HV_pre)./2)  nanstd((L2_c+L2_c)./2-(LA_pre+LV_pre)./2)]./sqrt(length(Pp)),'color',[0.4 0.5 0],'marker','x','linestyle',':','linewidth',6), hold on

% errorbar([1.7 3.7],[ nanmean((Hn_c+Hn_c)./2-(HA_pre+HV_pre)./2) nanmean((Ln_c+Ln_c)./2-(LA_pre+LV_pre)./2)  ],...
%     [nanstd((Hn_c+Hn_c)./2-(HA_pre+HV_pre)./2)  nanstd((Ln_c+Ln_c)./2-(LA_pre+LV_pre)./2)  ]./sqrt(length(Pp)),'color',[0.5 0.5 0.5],'marker','d','linewidth',1), hold on


errorbar([1.3 3.3]+5,[ nanmean((H1_e+H1_e)./2-(HA_pre+HV_pre)./2) nanmean((L1_e+L1_e)./2-(LA_pre+LV_pre)./2)  ],...
    [nanstd((H1_e+H1_e)./2-(HA_pre+HV_pre)./2)  nanstd((L1_e+L1_e)./2-(LA_pre+LV_pre)./2)  ]./sqrt(length(Pp)),'color',[1 0.5 0],'marker','d','linewidth',6), hold on
errorbar([1.5 3.5]+5,[ nanmean((H2_e+H2_e)./2-(HA_pre+HV_pre)./2)  nanmean((L2_e+L2_e)./2-(LA_pre+LV_pre)./2) ],...
    [nanstd((H2_e+H2_e)./2-(HA_pre+HV_pre)./2)  nanstd((L2_e+L2_e)./2-(LA_pre+LV_pre)./2)]./sqrt(length(Pp)),'color',[0.4 0.5 0],'marker','x','linestyle',':','linewidth',6), hold on
% errorbar([1.7 3.7]+5,[ nanmean((Hn_e+Hn_e)./2-(HA_pre+HV_pre)./2) nanmean((Ln_e+Ln_e)./2-(LA_pre+LV_pre)./2)  ],...
%     [nanstd((Hn_e+Hn_e)./2-(HA_pre+HV_pre)./2)  nanstd((Ln_e+Ln_e)./2-(LA_pre+LV_pre)./2)  ]./sqrt(length(Pp)),'color',[0.5 0.5 0.5],'marker','d','linewidth',1), hold on



% legend('after-High','after-Low','Location','SouthEast')
text(1, myYLim(2),'Trial n-1 Correct','FontSize',24)
text(6.5, myYLim(2),'Trial n-1 Error','FontSize',24)

ax = gca;
ax.XTickLabel = {'High','Low','High','Low'};
ax.YLim = myYLim;
ax.XTick = [ 1  3  6 8 ];
plot([0 9],[0 0],':k')
xlabel('Reward Current Trial')
ylabel('Accuracy')
[h p_hlce1]= ttest((H1_c-(HA_pre+HV_pre)./2)-(L1_c-(LA_pre+LV_pre)./2),(H2_c-(HA_pre+HV_pre)./2)-(L2_c-(LA_pre+LV_pre)./2))
[h p_hlce2]= ttest((H1_e-(HA_pre+HV_pre)./2)-(L1_e-(LA_pre+LV_pre)./2),(H2_e-(HA_pre+HV_pre)./2)-(L2_e-(LA_pre+LV_pre)./2))
[h p_hlce3]= ttest((H1_c-(HA_pre+HV_pre)./2)-(L1_c-(LA_pre+LV_pre)./2)-(H2_c-(HA_pre+HV_pre)./2)-(L2_c-(LA_pre+LV_pre)./2),(H1_e-(HA_pre+HV_pre)./2)-(L1_e-(LA_pre+LV_pre)./2)-(H2_e-(HA_pre+HV_pre)./2)-(L2_e-(LA_pre+LV_pre)./2))
[h p_hlce2]= ttest((H1_e-(HA_pre+HV_pre)./2)+(H2_e-(HA_pre+HV_pre)./2), (L1_e-(LA_pre+LV_pre)./2)+(L2_e-(LA_pre+LV_pre)./2))
[h p_hlce3]= ttest((H2_c-(HA_pre+HV_pre)./2)-(L2_c-(LA_pre+LV_pre)./2), (H2_e-(HA_pre+HV_pre)./2)-(L2_e-(LA_pre+LV_pre)./2))
[h p_hlce4]= ttest((Hn_c-(HA_pre+HV_pre)./2)-(Ln_c-(LA_pre+LV_pre)./2), (Hn_e-(HA_pre+HV_pre)./2)-(Ln_e-(LA_pre+LV_pre)./2))
%%
anovdata_rew_modality_names = {'HV1_c';'LV1_c';'HA1_c';'LA1_c';...
    'HV2_c';'LV2_c';'HA2_c';'LA2_c';...
    'HV1_e';'LV1_e';'HA1_e';'LA1_e';...
    'HV2_e';'LV2_e';'HA2_e';'LA2_e'};

% HV_pre= zeros(length(Pp),1)';
% LV_pre= zeros(length(Pp),1)';
% HA_pre= zeros(length(Pp),1)';
% LA_pre= zeros(length(Pp),1)';
% N_pre= zeros(length(Pp),1)';
%

anovdata_rew_modality = [HV1_c-HV_pre;LV1_c-LV_pre;HA1_c-HA_pre;LA1_c-LA_pre;...
    HV2_c-HV_pre;LV2_c-LV_pre;HA2_c-HA_pre;LA2_c-LA_pre;...
    HV1_e-HV_pre;LV1_e-LV_pre;HA1_e-HA_pre;LA1_e-LA_pre;...
    HV2_e-HV_pre;LV2_e-LV_pre;HA2_e-HA_pre;LA2_e-LA_pre]';

bimodncond=length(anovdata_rew_modality_names);
subfactor=cell(bimodncond,1);

for  i=1:bimodncond
    subfactor{i}= ['F' num2str(i)];
end

t = array2table([anovdata_rew_modality (1:length(Pp))'],'VariableNames',[subfactor;{'Sujs'}]);
%factorNames = {'modality','reward','latencybin','part'};
factorNames = {'modality','reward','posthl','cor_err'}; % for each part


within = table({'V';'V';'A';'A';'V';'V';'A';'A';'V';'V';'A';'A';'V';'V';'A';'A'},... %modality: visual, auditory, bimodal
    {'H';'L';'H';'L';'H';'L';'H';'L';'H';'L';'H';'L';'H';'L';'H';'L'},... %reward: high or low
    {'post_h';'post_h';'post_h';'post_h';'post_l';'post_l';'post_l';'post_l';'post_h';'post_h';'post_h';'post_h';'post_l';'post_l';'post_l';'post_l'},... %phase: pre or post
    {'cor';'cor';'cor';'cor';'cor';'cor';'cor';'cor';'err';'err';'err';'err';'err';'err';'err';'err'},... %phase: pre or post
    'VariableNames',factorNames);
%
% fit the repeated measures model
rm = fitrm(t,['F1-F' num2str(bimodncond) '~1'],'WithinDesign',within);

% run my repeated measures anova here
%[ranovatbl] = ranova(rm, 'WithinModel','modality*reward*latencybin*part')
[ranovatbl] = ranova(rm, 'WithinModel','modality*reward*posthl*cor_err')
%% Also getting Neutral condition
anovdata_rew_modality_names = {'HV1_c';'LV1_c';'HA1_c';'LA1_c';...
    'HV2_c';'LV2_c';'HA2_c';'LA2_c';...
    'HV1_e';'LV1_e';'HA1_e';'LA1_e';...
    'HV2_e';'LV2_e';'HA2_e';'LA2_e';...
    'HVn_c';'LVn_c';'HAc_c';'Ln_c';...
    'HVn_e';'LVn_e';'HAn_e';'LAn_e'};
anovdata_rew_modality = [HV1_c-HV_pre;LV1_c-LV_pre;HA1_c-HA_pre;LA1_c-LA_pre;...
    HV2_c-HV_pre;LV2_c-LV_pre;HA2_c-HA_pre;LA2_c-LA_pre;...
    HV1_e-HV_pre;LV1_e-LV_pre;HA1_e-HA_pre;LA1_e-LA_pre;...
    HV2_e-HV_pre;LV2_e-LV_pre;HA2_e-HA_pre;LA2_e-LA_pre;...
    HVn_c-HV_pre;LVn_c-LV_pre;HAn_c-HA_pre;LAn_c-LA_pre;...
    HVn_e-HV_pre;LVn_e-LV_pre;HAn_e-HA_pre;LAn_e-LA_pre]';

bimodncond=length(anovdata_rew_modality_names);
subfactor=cell(bimodncond,1);

for  i=1:bimodncond
    subfactor{i}= ['F' num2str(i)];
end

t = array2table([anovdata_rew_modality (1:length(Pp))'],'VariableNames',[subfactor;{'Sujs'}]);
%factorNames = {'modality','reward','latencybin','part'};
factorNames = {'modality','reward','posthl','cor_err'}; % for each part


within = table({'V';'V';'A';'A';'V';'V';'A';'A';'V';'V';'A';'A';'V';'V';'A';'A';'V';'V';'A';'A';'V';'V';'A';'A'},... %modality: visual, auditory, bimodal
    {'H';'L';'H';'L';'H';'L';'H';'L';'H';'L';'H';'L';'H';'L';'H';'L';'H';'L';'H';'L';'H';'L';'H';'L'},... %reward: high or low
    {'post_h';'post_h';'post_h';'post_h';'post_l';'post_l';'post_l';'post_l';'post_h';'post_h';'post_h';'post_h';'post_l';'post_l';'post_l';'post_l';...
    'post_n';'post_n';'post_n';'post_n';'post_n';'post_n';'post_n';'post_n'},... %phase: pre or post
    {'cor';'cor';'cor';'cor';'cor';'cor';'cor';'cor';'err';'err';'err';'err';'err';'err';'err';'err';'cor';'cor';'cor';'cor';'err';'err';'err';'err'},... %phase: pre or post
    'VariableNames',factorNames);
%
% fit the repeated measures model
rm = fitrm(t,['F1-F' num2str(bimodncond) '~1'],'WithinDesign',within);

% run my repeated measures anova here
%[ranovatbl] = ranova(rm, 'WithinModel','modality*reward*latencybin*part')
[ranovatbl] = ranova(rm, 'WithinModel','modality*reward*posthl*cor_err')



%%
% pre=load('pre_beh5_jessica.mat','s','Pp');
% HV= s.subBCvH(Pp)';
% LV= s.subBCvL(Pp)';
% HA= s.subSPvH(Pp)';
% LA= s.subSPvL(Pp)';
% N= s.subNeut(Pp)';
% HV_pre= pre.s.subBCvH(Pp)';
% LV_pre= pre.s.subBCvL(Pp)';
% HA_pre= pre.s.subSPvH(Pp)';
% LA_pre= pre.s.subSPvL(Pp)';
% N_pre= pre.s.subNeut(Pp)';
% anovdata_rew_modality_names = {'HV';'LV';'HA';'LA'};
% anovdata_rew_modality = [HV;LV;HA;LA]';
%
% bimodncond=length(anovdata_rew_modality_names);
% subfactor=cell(bimodncond,1);
%
% for  i=1:bimodncond
%     subfactor{i}= ['F' num2str(i)];
% end
% accufac = zeros(length(Pp),1);
% for p =1:length(Pp)
%     if s.subDPri(Pp(p))>median(s.subDPri(Pp))
%         accufac(p)=1;
%     end
% end
%
% t = array2table([anovdata_rew_modality s.subRewS(Pp,1) s.subRewV(Pp,1) (accufac)],...
%     'VariableNames',[subfactor; {'SoundP'} ;{'BoxColor'};{'Accuracy'}]);
% t.SoundP = categorical (t.SoundP);
% t.BoxColor = categorical (t.BoxColor);
% t.Accuracy = categorical (t.Accuracy);
%
% factorNames = {'modality','reward'}; % for each part
%
%
% within = table({'V';'V';'A';'A'},... %modality: visual, auditory, bimodal
%     {'H';'L';'H';'L'},... %reward: high or low
%     'VariableNames',factorNames);
%
%
% % % fit the repeated measures model
% % rm = fitrm(t,['F1-F' num2str(bimodncond) '~1'],'WithinDesign',within);
%
% % fit the repeated measures model
% rm = fitrm(t,['F1-F' num2str(bimodncond) '~Accuracy'],'WithinDesign',within);
%
% % run my repeated measures anova here
% %[ranovatbl] = ranova(rm, 'WithinModel','modality*reward*latencybin*part')
% [ranovatbl] = ranova(rm, 'WithinModel','modality*reward')
%%
% pre=load('pre_beh5_jessica.mat','s','Pp');
% HV= s.subBCvH(Pp)';
% LV= s.subBCvL(Pp)';
% HA= s.subSPvH(Pp)';
% LA= s.subSPvL(Pp)';
% N= s.subNeut(Pp)';
% HV_pre= pre.s.subBCvH(Pp)';
% LV_pre= pre.s.subBCvL(Pp)';
% HA_pre= pre.s.subSPvH(Pp)';
% LA_pre= pre.s.subSPvL(Pp)';
% N_pre= pre.s.subNeut(Pp)';
% anovdata_rew_modality_names = {'HV';'LV';'HA';'LA';'HV_pre';'LV_pre';'HA_pre';'LA_pre'};
% anovdata_rew_modality = [HV;LV;HA;LA;HV_pre;LV_pre;HA_pre;LA_pre]';
%
% bimodncond=length(anovdata_rew_modality_names);
% subfactor=cell(bimodncond,1);
%
% for  i=1:bimodncond
%     subfactor{i}= ['F' num2str(i)];
% end
% accufac = zeros(length(Pp),1);
% for p =1:length(Pp)
%     if pre.s.subAccu(Pp(p))>median(pre.s.subAccu(Pp))
%         accufac(p)=1;
%     end
% end
%
% t = array2table([anovdata_rew_modality s.subRewS(Pp,1) s.subRewV(Pp,1) (s.subDPri(Pp,1))],...
%     'VariableNames',[subfactor; {'SoundP'} ;{'BoxColor'};{'Accuracy'}]);
% t.SoundP = categorical (t.SoundP);
% t.BoxColor = categorical (t.BoxColor);
% %t.Accuracy = categorical (t.Accuracy);
%
% factorNames = {'modality','reward','pre_post'}; % for each part
%
%
% within = table({'V';'V';'A';'A';'V';'V';'A';'A'},... %modality: visual, auditory, bimodal
%     {'H';'L';'H';'L';'H';'L';'H';'L'},... %reward: high or low
%      {'Post';'Post';'Post';'Post';'Pre';'Pre';'Pre';'Pre'},... %phase: pre or post
%     'VariableNames',factorNames);
%
%
% % % fit the repeated measures model
% % rm = fitrm(t,['F1-F' num2str(bimodncond) '~1'],'WithinDesign',within);
%
% % fit the repeated measures model
% rm = fitrm(t,['F1-F' num2str(bimodncond) '~Accuracy'],'WithinDesign',within);
%
% % run my repeated measures anova here
% %[ranovatbl] = ranova(rm, 'WithinModel','modality*reward*latencybin*part')
% [ranovatbl] = ranova(rm, 'WithinModel','modality*reward*pre_post')
%%
%% Figure 2 Sound
% %% Figure 2 Sound
% pre=load('pre_beh5_jessica.mat','s','Pp');
% HV= s.subBCvH(Pp)';
% LV= s.subBCvL(Pp)';
% HA= s.subSPvH(Pp)';
% LA= s.subSPvL(Pp)';
% N= s.subNeut(Pp)';
% HV_pre= pre.s.subBCvH(Pp)';
% LV_pre= pre.s.subBCvL(Pp)';
% HA_pre= pre.s.subSPvH(Pp)';
% LA_pre= pre.s.subSPvL(Pp)';
% N_pre= pre.s.subNeut(Pp)';
% figure(1700), subplot(2,1,1), hold on;
% bar([1 3 4 6 7],[nanmean(N-N_pre) nanmean(HV-HV_pre) nanmean(LV-LV_pre) nanmean(HA-HA_pre)...
%     nanmean(LA-LA_pre)],'w'), axis square, hold on
% errorbar([1 3 4 6 7],[nanmean(N-N_pre) nanmean(HV-HV_pre) nanmean(LV-LV_pre) nanmean(HA-HA_pre)...
%     nanmean(LA-LA_pre)],...
%     [nanstd(N-N_pre) nanstd(HV-HV_pre) nanstd(LV-LV_pre) nanstd(HA-HA_pre)...
%     nanstd(LA-LA_pre)]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none')
% ax = gca;
% ax.XTickLabel = {'Neut','VisHighVal','VisLowVal','AudioHighVal','AudioLowVal'};
% % if whichone==1
% %     ax.YLim = scalesel1;
% %     ylabel('DPrime' , 'FontSize', 14);
% % elseif whichone==2
% %     ax.YLim = scalesel2;
% %     ylabel('Accuracy' , 'FontSize', 14);
% % elseif whichone==3
% %     ax.YLim = scalesel3;
% %     ylabel('RT' , 'FontSize', 14);
% % elseif whichone==5
% %     ax.YLim = scalesel5;
% %     ylabel('Efficiency' , 'FontSize', 14);
% % elseif whichone==6
% %     ax.YLim = scalesel6;
% %     ylabel('Key (1-UP, 2-Down)' , 'FontSize', 14);
% % end
% ax = gca;
% % ax.YLim = [1.4 1.6];
% title('Effect of Reward value and modality', 'FontSize', 16);
% [h p1]= ttest(HV-HV_pre,LV-LV_pre);
% [h p2]= ttest(HA-HA_pre,LA-LA_pre);
% [h p3]= ttest((HV+HA)-(HV_pre+HA_pre),(LV+LA)-(LV_pre+LA_pre))
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', visual P= ' num2str(round(p1,3)) ', auditory P= ' num2str(round(p2,3))]);
% % print ('-f2' ,'-dpsc2',[rootdir_fig 'valsounds.ps'])
% %%
% pre=load('pre_beh5_jessica.mat','s','Pp');
% Accu = s.subDPri(Pp);
% pre_Accu= pre.s.subDPri(Pp);
%
% figure(1700), subplot(2,1,2), hold on;
% bar([1 3 ],[nanmean(Accu) nanmean(pre_Accu) ],'w'), axis square, hold on
% errorbar([1 3 ],[nanmean(Accu) nanmean(pre_Accu) ],...
%     [nanstd(Accu) nanstd(pre_Accu) ]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none')
% ax = gca;
% ax.XTick = [0.5 1];
% ax.XTickLabel = {'Accu','Pre-Accu'};
%
% ax = gca;
% ax.YLim = [0.5 1 ];
% title('Learning', 'FontSize', 16);
% [h p1]= ttest(Accu,pre_Accu);
% [h p2]= ttest(HA-HA_pre,LA-LA_pre);
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', Learning P= ' num2str(round(p1,3)) ]);
% % print ('-f2' ,'-dpsc2',[rootdir_fig 'valsounds.ps'])
%%
pre=load('pre_beh5_jessica.mat','s','Pp');
post1=load('post_beh5_jessica_1sthalf.mat','s','Pp');
post2=load('post_beh5_jessica_2ndhalf.mat','s','Pp');

HV1= post1.s.subBCvH(Pp)';
LV1= post1.s.subBCvL(Pp)';
HA1= post1.s.subSPvH(Pp)';
LA1= post1.s.subSPvL(Pp)';
N1= post1.s.subNeut(Pp)';

HV2= post2.s.subBCvH(Pp)';
LV2= post2.s.subBCvL(Pp)';
HA2= post2.s.subSPvH(Pp)';
LA2= post2.s.subSPvL(Pp)';
N2= post2.s.subNeut(Pp)';


HV_pre= pre.s.subBCvH(Pp)';
LV_pre= pre.s.subBCvL(Pp)';
HA_pre= pre.s.subSPvH(Pp)';
LA_pre= pre.s.subSPvL(Pp)';
N_pre= pre.s.subNeut(Pp)';

xloci =[1 11 20 31 40];
figure(1710),  hold on;
bar([xloci xloci+xloci(end)+12],[nanmean(N1-N_pre) nanmean(HV1-HV_pre) nanmean(LV1-LV_pre) nanmean(HA1-HA_pre)...
    nanmean(LA1-LA_pre) nanmean(N2-N_pre) nanmean(HV2-HV_pre) nanmean(LV2-LV_pre) nanmean(HA2-HA_pre)...
    nanmean(LA2-LA_pre)],1,'w'), axis square, hold on
errorbar([xloci xloci+xloci(end)+12],[nanmean(N1-N_pre) nanmean(HV1-HV_pre) nanmean(LV1-LV_pre) nanmean(HA1-HA_pre)...
    nanmean(LA1-LA_pre) nanmean(N2-N_pre) nanmean(HV2-HV_pre) nanmean(LV2-LV_pre) nanmean(HA2-HA_pre)...
    nanmean(LA2-LA_pre)],...
    ([nanstd(N1-N_pre) nanstd(HV1-HV_pre) nanstd(LV1-LV_pre) nanstd(HA1-HA_pre)...
    nanstd(LA1-LA_pre) nanstd(N2-N_pre) nanstd(HV2-HV_pre) nanstd(LV2-LV_pre) nanstd(HA2-HA_pre)...
    nanstd(LA2-LA_pre)])./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none')
ax = gca;
ax.XTickLabel = {'Neut','VisHighVal','VisLowVal','AudioHighVal','AudioLowVal'};
% if whichone==1
%     ax.YLim = scalesel1;
%     ylabel('DPrime' , 'FontSize', 14);
% elseif whichone==2
%     ax.YLim = scalesel2;
%     ylabel('Accuracy' , 'FontSize', 14);
% elseif whichone==3
%     ax.YLim = scalesel3;
%     ylabel('RT' , 'FontSize', 14);
% elseif whichone==5
%     ax.YLim = scalesel5;
%     ylabel('Efficiency' , 'FontSize', 14);
% elseif whichone==6
%     ax.YLim = scalesel6;
%     ylabel('Key (1-UP, 2-Down)' , 'FontSize', 14);
% end
ax = gca;
%ax.YLim = [-0.05 .1];
title('Effect of Reward value and modality', 'FontSize', 16);
[h p1]= ttest(HV1-HV_pre,LV1-LV_pre);
[h p2]= ttest(HA1-HA_pre,LA1-LA_pre);
[h p3]= ttest((HV1+HA1)-(HV_pre+HA_pre),(LV1+LA1)-(LV_pre+LA_pre))

[h p4]= ttest(HV2-HV_pre,LV2-LV_pre);
[h p5]= ttest(HA2-HA_pre,LA2-LA_pre);
[h p6]= ttest((HV2+HA2)-(HV_pre+HA_pre),(LV2+LA2)-(LV_pre+LA_pre))
title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
    ', N= ' num2str(length(Pp)) ', visual P= ' num2str(round(p1,3)) ', auditory P= ' num2str(round(p2,3))  ', visual2 P= ' num2str(round(p4,3)) ', auditory2 P= ' num2str(round(p5,3))]);
% print ('-f2' ,'-dpsc2',[rootdir_fig 'valsounds.ps'])
%%
% figure(1720),errorbar([1 2 ],[nanmean(N1-N_pre)  nanmean(N2-N_pre)],...
%     [nanstd(N1-N_pre) nanstd(N2-N_pre) ]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle',':'), hold on
% errorbar([1.01 2.01 ],[nanmean(HV1-HV_pre)  nanmean(HV2-HV_pre)],...
%     [nanstd(HV1-HV_pre) nanstd(HV2-HV_pre) ]./sqrt(length(Pp)),'r','linewidth',3, 'linestyle',':')
% errorbar([1.02 2.02 ],[nanmean(LV1-LV_pre)  nanmean(LV2-LV_pre)],...
%     [nanstd(LV1-LV_pre) nanstd(LV2-LV_pre) ]./sqrt(length(Pp)),'b','linewidth',3, 'linestyle',':')
% errorbar([1.03 2.03 ],[nanmean(HA1-HA_pre)  nanmean(HA2-HA_pre)],...
%     [nanstd(HA1-HA_pre) nanstd(HA2-HA_pre) ]./sqrt(length(Pp)),'m','linewidth',3, 'linestyle',':')
% errorbar([1.04 2.04 ],[nanmean(LA1-LA_pre)  nanmean(LA2-LA_pre)],...
%     [nanstd(LA1-LA_pre) nanstd(LA2-LA_pre) ]./sqrt(length(Pp)),'c','linewidth',3, 'linestyle',':')
% %%
% figure(1740),errorbar([1 2 ],[nanmean((HV1-LV1))  nanmean(HV2-LV2)],...
%     [nanstd(HV1-LV1) nanstd(HV2-LV2) ]./sqrt(length(Pp)),'r','linewidth',3, 'linestyle',':'), hold on
% errorbar([1.01 2.01 ],[nanmean(HA1-LA1)  nanmean(HA2-LA2)],...
%     [nanstd(HA1-LA1) nanstd(HA2-LA2) ]./sqrt(length(Pp)),'b','linewidth',3, 'linestyle',':')
% %%
% figure(1750),errorbar([1 2 ],[nanmean((HV1-LV1)-(HV_pre-LV_pre))  nanmean((HV2-LV2)-(HV_pre-LV_pre))],...
%     [nanstd((HV1-LV1)-(HV_pre-LV_pre)) nanstd((HV2-LV2)-(HV_pre-LV_pre)) ]./sqrt(length(Pp)),'r','linewidth',3, 'linestyle',':'), hold on
% errorbar([1.1 2.1 ],[nanmean((HA1-LA1)-(HA_pre-LA_pre))  nanmean((HA2-LA2)-(HA_pre-LA_pre))],...
%     [nanstd((HA1-LA1)-(HA_pre-LA_pre)) nanstd((HA2-LA2)-(HA_pre-LA_pre)) ]./sqrt(length(Pp)),'b','linewidth',3, 'linestyle',':'), hold on
%%
pre=load('pre_beh5_jessica.mat','s','Pp');

HV = s.bin_subBCvH(Pp,:);
LV = s.bin_subBCvL(Pp,:);
HA = s.bin_subSPvH(Pp,:);
LA = s.bin_subSPvL(Pp,:);
N = s.bin_subNeut(Pp,:);

HV_pre= pre.s.subBCvH(Pp)';
LV_pre= pre.s.subBCvL(Pp)';
HA_pre= pre.s.subSPvH(Pp)';
LA_pre= pre.s.subSPvL(Pp)';
N_pre= pre.s.subNeut(Pp)';
%%%%% if we want to normalize
% for i=1:length(Pp)
%     maxi(i) = max([max(HV(i,:));max(LV(i,:));max(HA(i,:));max(LA(i,:));max(N(i,:));max(HV_pre(i));max(LV_pre(i));max(HA_pre(i));max(LA_pre(i));max(N_pre(i))]);
%     HV(i,:) = HV(i,:)./maxi(i);
%     LV(i,:) = LV(i,:)./maxi(i);
%     HA(i,:) = HA(i,:)./maxi(i);
%     LA(i,:) = LA(i,:)./maxi(i);
%     N(i,:)  = N(i,:)./maxi(i);
%
%     HV_pre(i) = HV_pre(i)'./maxi(i);
%     LV_pre(i) = LV_pre(i)'./maxi(i);
%     HA_pre(i) = HA_pre(i)'./maxi(i);
%     LA_pre(i) = LA_pre(i)'./maxi(i);
%     N_pre(i)  = N_pre(i)'./maxi(i);
% end

figure (888),axis square, hold on

errorbar([1:nbins-1], [mean(N-repmat(N_pre',1,nbins-1)) ],[std(N-repmat(N_pre',1,nbins-1)) ]./sqrt(length(Pp)),'k','linewidth',3), hold on
errorbar([1:nbins-1], [mean(HV-repmat(HV_pre',1,nbins-1)) ],[std(HV-repmat(HV_pre',1,nbins-1)) ]./sqrt(length(Pp)),'r','linewidth',3), hold on
errorbar([1:nbins-1], [mean(LV-repmat(LV_pre',1,nbins-1)) ],[std(LV-repmat(LV_pre',1,nbins-1)) ]./sqrt(length(Pp)),':m','linewidth',3), hold on
errorbar([1:nbins-1], [mean(HA-repmat(HA_pre',1,nbins-1)) ],[std(HA-repmat(HA_pre',1,nbins-1)) ]./sqrt(length(Pp)),'b','linewidth',3), hold on
errorbar([1:nbins-1], [mean(LA-repmat(LA_pre',1,nbins-1)) ],[std(LA-repmat(LA_pre',1,nbins-1)) ]./sqrt(length(Pp)),':c','linewidth',3), hold on

ax = gca;
ax.YLim = [-0.05 .15];

%
figure (889),axis square, hold on
errorbar([1:nbins-1], [mean((HV-repmat(HV_pre',1,nbins-1))-(LV-repmat(LV_pre',1,nbins-1))) ],...
    [std((HV-repmat(HV_pre',1,nbins-1))-(LV-repmat(LV_pre',1,nbins-1))) ]./sqrt(length(Pp)),'color',[0.1 0.5 0.5],'linewidth',3,'marker','x','markersize',8), hold on
errorbar([1:nbins-1]+0.05, [mean((HA-repmat(HA_pre',1,nbins-1))-(LA-repmat(LA_pre',1,nbins-1))) ],...
    [std((HA-repmat(HA_pre',1,nbins-1))-(LA-repmat(LA_pre',1,nbins-1))) ]./sqrt(length(Pp)),'color',[0.5 0.1 0.5],'linewidth',3,'marker','o','markersize',8), hold on
plot([ 0 4],[0 0],':k')
ax = gca;
ax.YLim = [-0.1 .1];

% figure (890), hold on
% errorbar([1:nbins-1], [mean((HV)-(LV)) ],...
%     [std((HV)-(LV)) ]./sqrt(length(Pp)),'r','linewidth',3), hold on
% errorbar([1:nbins-1]+0.2, [mean((HA)-(LA)) ],...
%     [std((HA)-(LA)) ]./sqrt(length(Pp)),'b','linewidth',3), hold on


% figure (900), hold on
% errorbar([1:nbins-1], [mean((HV-repmat(HV_pre',1,nbins-1))-(N-repmat(N_pre',1,nbins-1))) ],...
%     [std((HV-repmat(HV_pre',1,nbins-1))-(N-repmat(N_pre',1,nbins-1))) ]./sqrt(length(Pp)),'r','linewidth',3), hold on
% errorbar([1:nbins-1], [mean((LV-repmat(LV_pre',1,nbins-1))-(N-repmat(N_pre',1,nbins-1))) ],...
%     [std((LV-repmat(LV_pre',1,nbins-1))-(N-repmat(N_pre',1,nbins-1))) ]./sqrt(length(Pp)),':m','linewidth',3), hold on
[h p11]=ttest((HA-repmat(HA_pre',1,nbins-1)),(LA-repmat(LA_pre',1,nbins-1)))
[h p12]=ttest((HV-repmat(HV_pre',1,nbins-1)),(LV-repmat(LV_pre',1,nbins-1)))
[h p13]= ttest((HA-repmat(HA_pre',1,nbins-1))+(HV-repmat(HV_pre',1,nbins-1)),(LA-repmat(LA_pre',1,nbins-1))+(LV-repmat(LV_pre',1,nbins-1))) %%% for interaction
title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
    ', N= ' num2str(length(Pp)) ', Learning-diff-highlow P= ' num2str(round(p13,3)) ]);


%% Looking at the n-1 and n trials dependent on the reward

pre=load('pre_beh5_jessica.mat','s','Pp');

HV = s.bin_subBCvH(Pp,:);
LV = s.bin_subBCvL(Pp,:);
HA = s.bin_subSPvH(Pp,:);
LA = s.bin_subSPvL(Pp,:);
N = s.bin_subNeut(Pp,:);

HV_pre= pre.s.subBCvH(Pp)';
LV_pre= pre.s.subBCvL(Pp)';
HA_pre= pre.s.subSPvH(Pp)';
LA_pre= pre.s.subSPvL(Pp)';
N_pre= pre.s.subNeut(Pp)';

% HV_pre= zeros(length(Pp),1)';
% LV_pre= zeros(length(Pp),1)';
% HA_pre= zeros(length(Pp),1)';
% LA_pre= zeros(length(Pp),1)';
% N_pre= zeros(length(Pp),1)';



% figure(17320), axis square, hold on;
% boxplot([((HV+HA)./2-(repmat(HV_pre',1,nbins-1)+repmat(HA_pre',1,nbins-1))./2)  ((LV+LA)./2-(repmat(LV_pre',1,nbins-1)+repmat(LA_pre',1,nbins-1))./2) (N-repmat(N_pre',1,nbins-1))])


figure(17330), axis square, hold on;
% bar([1 ],[ nanmean((HV+HA)./2-(HV_pre+HA_pre)./2)  ],'r'), hold on
%
% bar([3],[  nanmean((LV+LA)./2-(LV_pre+LA_pre)./2) ],'b')
%
% bar([5 ],[  nanmean(N-N_pre)],'k')

errorbar([1.1 2.1], nanmean(HV-repmat(HV_pre',1,nbins-1)),nanstd(HV-repmat(HV_pre',1,nbins-1))./sqrt(length(Pp)),':c','linewidth',3, 'marker','*'), hold on
errorbar([1.2 2.2], nanmean(HA-repmat(HA_pre',1,nbins-1)),nanstd(HA-repmat(HA_pre',1,nbins-1))./sqrt(length(Pp)),':m','linewidth',3, 'marker','*'), hold on

errorbar([1 2],[ nanmean((HV+HA)./2-(repmat(HV_pre',1,nbins-1)+repmat(HA_pre',1,nbins-1))./2)  ],...
    [nanstd((HV+HA)./2-(repmat(HV_pre',1,nbins-1)+repmat(HA_pre',1,nbins-1))./2)  ]./sqrt(length(Pp)),'r','linewidth',3, 'marker','d'), hold on

errorbar([3.1 4.1], nanmean(LV-repmat(LV_pre',1,nbins-1)),nanstd(LV-repmat(LV_pre',1,nbins-1))./sqrt(length(Pp)),':c','linewidth',3, 'marker','*'), hold on
errorbar([3.2 4.2], nanmean(LA-repmat(LA_pre',1,nbins-1)),nanstd(LA-repmat(LA_pre',1,nbins-1))./sqrt(length(Pp)),':m','linewidth',3, 'marker','*'), hold on
errorbar([3 4],[  nanmean((LV+LA)./2-(repmat(LV_pre',1,nbins-1)+repmat(LA_pre',1,nbins-1))./2) ],...
    [ nanstd((LV+LA)./2-(repmat(LV_pre',1,nbins-1)+repmat(LA_pre',1,nbins-1))./2)  ]./sqrt(length(Pp)),'b','linewidth',3, 'marker','d')

errorbar([5 6 ],[  nanmean(N-repmat(N_pre',1,nbins-1))],...
    [ nanstd(N-repmat(N_pre',1,nbins-1))  ]./sqrt(length(Pp)),'k','linewidth',3, 'marker','d')

ax = gca;
ax.XTickLabel = {'High','Low','Neutral'};
ax.XTick = [ 1  3  5 ];
ax.YLim = [-.05 0.1];
plot([0 6],[0 0],':k')
[h p_hh]= ttest((HV+HA)./2,(repmat(HV_pre',1,nbins-1)+repmat(HA_pre',1,nbins-1))./2)  %%%%% this is switch versus no switch
[h p_ll]= ttest((LV+LA)./2,(repmat(LV_pre',1,nbins-1)+repmat(LA_pre',1,nbins-1))./2)  %%%%% this is switch versus no switch
[h p_n]= ttest(N,repmat(N_pre',1,nbins-1))  %%%%% this is switch versus no switch
[h p_hl]= ttest(((HV+HA)./2-(repmat(HV_pre',1,nbins-1)+repmat(HA_pre',1,nbins-1))./2),((LV+LA)./2-(repmat(LV_pre',1,nbins-1)+repmat(LA_pre',1,nbins-1))./2 ))
%[h p_hl]= ttest((HV+HA)./2-(HV_pre+HA_pre)./2,N-N_pre)


title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
    ', N= ' num2str(length(Pp)) ', Learning_H P= ' num2str(round(p_hh,3)) ', Learning_L P= ' num2str(round(p_ll,3))...
    ', Learning_N P= ' num2str(round(p_n,3)) ', Difference HIgh and Low P= ' num2str(round(p_hl,3))]);


%%
pre=load('pre_beh5_jessica.mat','s','Pp');
post1H=load('post_beh5_jessica_afterHigh_nminus1_correct.mat','s','Pp');
post2H=load('post_beh5_jessica_afterHigh2_correct.mat','s','Pp');

post1L=load('post_beh5_jessica_afterLow_nminus1_correct.mat','s','Pp');
post2L=load('post_beh5_jessica_afterLow2_correct.mat','s','Pp');


% pre=load('pre_beh5_jessica.mat','s','Pp');
% post1H=load('post_beh5_jessica_afterHigh_nminus1.mat','s','Pp');
% post2H=load('post_beh5_jessica_afterHigh2.mat','s','Pp');
%
% post1L=load('post_beh5_jessica_afterLow_nminus1.mat','s','Pp');
% post2L=load('post_beh5_jessica_afterLow2.mat','s','Pp');


n1h= (post1H.s.subBCvH(Pp)'+post1H.s.subBCvL(Pp)'+post1H.s.subSPvH(Pp)'+ post1H.s.subSPvL(Pp)'+post1H.s.subNeut(Pp)')/5;

n1l= (post1L.s.subBCvH(Pp)'+post1L.s.subBCvL(Pp)'+post1L.s.subSPvH(Pp)'+ post1L.s.subSPvL(Pp)'+post1L.s.subNeut(Pp)')/5;


HV1= post2H.s.subBCvH(Pp)';
LV1= post2H.s.subBCvL(Pp)';
HA1= post2H.s.subSPvH(Pp)';
LA1= post2H.s.subSPvL(Pp)';
N1= post2H.s.subNeut(Pp)';

HV2= post2L.s.subBCvH(Pp)';
LV2= post2L.s.subBCvL(Pp)';
HA2= post2L.s.subSPvH(Pp)';
LA2= post2L.s.subSPvL(Pp)';
N2= post2L.s.subNeut(Pp)';

% %%%%% to test all trials
% HV1= s.subBCvH(Pp)';
% LV1= s.subBCvL(Pp)';
% HA1= s.subSPvH(Pp)';
% LA1= s.subSPvL(Pp)';
% N1= s.subNeut(Pp)';
%
% HV2= s.subBCvH(Pp)';
% LV2= s.subBCvL(Pp)';
% HA2= s.subSPvH(Pp)';
% LA2= s.subSPvL(Pp)';
% N2=  s.subNeut(Pp)';
% %%%to test all trials

HV_pre= pre.s.subBCvH(Pp)';
LV_pre= pre.s.subBCvL(Pp)';
HA_pre= pre.s.subSPvH(Pp)';
LA_pre= pre.s.subSPvL(Pp)';
N_pre= pre.s.subNeut(Pp)';

% HV_pre= zeros(length(Pp),1)';
% LV_pre= zeros(length(Pp),1)';
% HA_pre= zeros(length(Pp),1)';
% LA_pre= zeros(length(Pp),1)';
% N_pre= zeros(length(Pp),1)';


% figure(1810), subplot(2,1,1), hold on;
% bar([1 3 4 6 7 10 13 14 16 17],[nanmean(n1h) nanmean(HV1) nanmean(LV1) nanmean(HA1)...
%     nanmean(LA1) nanmean(n1l) nanmean(HV2) nanmean(LV2) nanmean(HA2)...
%     nanmean(LA2)],1,'w'), axis square, hold on
% errorbar([1 3 4 6 7 10 13 14 16 17],[nanmean(n1h) nanmean(HV1) nanmean(LV1) nanmean(HA1)...
%     nanmean(LA1) nanmean(n1l) nanmean(HV2) nanmean(LV2) nanmean(HA2)...
%     nanmean(LA2)],...
%     [nanstd(n1h) nanstd(HV1) nanstd(LV1) nanstd(HA1)...
%     nanstd(LA1) nanstd(n1l) nanstd(HV2) nanstd(LV2) nanstd(HA2)...
%     nanstd(LA2)]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none')
% ax = gca;
% ax.XTickLabel = {'Neut','VisHighVal','VisLowVal','AudioHighVal','AudioLowVal'};
% ax = gca;
% ax.YLim = [.6 1];
%
% figure(1820),  hold on;
% bar([1 3 6 8 11 13 16 18 21 23 ],[nanmean(n1h) nanmean(n1l) nanmean(HV1) nanmean(HV2) nanmean(LV2) nanmean(LV1) nanmean(HA1)...
%      nanmean(HA2) nanmean(LA2)  nanmean(LA1) ],'w'), axis square, hold on
% errorbar([1 3 6 8 11 13 16 18 21 23  ],[nanmean(n1h) nanmean(n1l) nanmean(HV1) nanmean(HV2) nanmean(LV2) nanmean(LV1) nanmean(HA1)...
%      nanmean(HA2) nanmean(LA2)  nanmean(LA1) ],...
%     [nanstd(n1h) nanstd(n1l) nanstd(HV1) nanstd(HV2) nanstd(LV2) nanstd(LV1) nanstd(HA1)...
%      nanstd(HA2) nanstd(LA2)  nanstd(LA1) ]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none')
% ax = gca;
% ax.XTick = [1 3 6 8 11 13 16 18 21 23  ];
% ax.XTickLabel = {'n-1h','n-1l','VHNS','VHS','VLNS','VLS','AHNS','AHS','ALNS','ALS'};
% ax.YLim = [.6 1];
%
% % title('Effect of Reward value and modality', 'FontSize', 16);
% % [h p1]= ttest(HV1,LV1);
% % [h p2]= ttest(HA1,LA1);
% % [h p3]= ttest((HV1+HA1)-(HV_pre+HA_pre),(LV1+LA1)-(LV_pre+LA_pre))
% %
% % [h p4]= ttest(HV2,LV2);
% % [h p5]= ttest(HA2,LA2);
% % [h p6]= ttest((HV2+HA2)-(HV_pre+HA_pre),(LV2+LA2)-(LV_pre+LA_pre))
% % title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
% %     ', N= ' num2str(length(Pp)) ', visual P= ' num2str(round(p1,3)) ', auditory P= ' num2str(round(p2,3))  ', visual2 P= ' num2str(round(p4,3)) ', auditory2 P= ' num2str(round(p5,3))]);
% % % print ('-f2' ,'-dpsc2',[rootdir_fig 'valsounds.ps'])
%
%
% figure(1821), subplot(2,1,1), hold on;
% bar([1:2:8],[nanmean((HV1+HA1)/2) nanmean((LV1+LA1)/2) nanmean((HV2+HA2)/2) nanmean((LV2+LA2)/2) ],'w'), axis square, hold on
% errorbar([1:2:8],[nanmean((HV1+HA1)/2) nanmean((LV1+LA1)/2) nanmean((HV2+HA2)/2) nanmean((LV2+LA2)/2) ],...
%     [nanstd((HV1+HA1)/2) nanstd((LV1+LA1)/2) nanstd((HV2+HA2)/2) nanstd((LV2+LA2)/2) ]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none')
% ax = gca;
% ax.XTick = [1:2:8 ];
% ax.XTickLabel = {'After H-H','After H-L','After L-H','After L-L'};
% ax.YLim = [.6 1];
% [h p1]= ttest((HV1+HA1)/2,(LV1+LA1)/2)
% [h p2]= ttest((HV2+HA2)/2,(LV2+LA2)/2)
% [h p3]= ttest((HV1+HA1)/2-(LV1+LA1)/2,(HV2+HA2)/2-(LV2+LA2)/2) %%%%5 interaction
% [h p4]= ttest((HV1+HA1+HV2+HA2)/4,(LV1+LA1+LV2+LA2)/4) %%%%all High,all low
%
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', After High P= ' num2str(round(p1,3)) ', After Low P= ' num2str(round(p2,3)) ', interaction P= ' num2str(round(p3,3)) ]);
%
% [h p10]= ttest((HA1-LA1)/2,(HA2-LA2)/2) %%%%% this is the interaction
% [h p11]= ttest((HV1-LV1)/2,(HV2-LV2)/2)
%
% [h p12]= ttest((HA1+HV1)/2,(HA2+HV2)/2)  %%%%% this is switch versus no switch
%
% [h p13]= ttest((LA1+LV1)/2,(LA2+LV2)/2)  %%%%% this is switch versus no switch
%
% [h p120]= ttest((HA1+HV1)/2-n1h,(HA2+HV2)/2-n1l)  %%%%% this is switch versus no switch
%
% [h p130]= ttest((LA1+LV1)/2-n1h,(LA2+LV2)/2-n1l)  %%%%% this is switch versus no switch
%
% [h p140]= ttest((LA1+LV1)/2-n1l,(HA1+HV1)/2-n1h)   %%%%% high no switch, low switch
%
% [h p150]= ttest((LA1+LV1)/2-n1l,(HA2+HV2)/2-n1l)   %%%%% high  switch, low switch
%
[h p1200]= ttest((HV1)-n1h,(HV2)-n1l)  %%%%% this is switch versus no switch

[h p1201]= ttest((HA1)-n1h,(HA2)-n1l)  %%%%% this is switch versus no switch

[h p1300]= ttest((LV1)-n1l,(LV2)-n1h)  %%%%% this is switch versus no switch
[h p1301]= ttest((LA1)-n1l,(LA2)-n1h)  %%%%% this is switch versus no switch
%
%
% [h p140]= ttest(((HA1+HV1)/2-n1h)-((HA2+HV2)/2-n1l),((LV2+LA2)/2-n1l)-((LA1+LV1)/2-n1h))  %%%%% this is interaction that the cost of switch is different for high and low reward
%
% figure(1830), axis square, hold on;
% errorbar([1 2],[ nanmean((HV1+HA1)./2) nanmean((HV2+HA2)./2) ],...
%     [nanstd((HV1+HA1)./2) nanstd((HV2+HA2)./2) ]./sqrt(length(Pp)),'r','linewidth',3), hold on
%
% errorbar([4 5],[  nanmean((LV2+LA2)./2) nanmean((LV1+LA1)./2)],...
%     [ nanstd((LV2+LA2)./2)  nanstd((LV1+LA1)./2)]./sqrt(length(Pp)),'b','linewidth',3)
%
% ax = gca;
% ax.XTickLabel = {'No-Switch-H','Switch-H','No-Switch-L','Switch-L'};
% ax.XTick = [ 1 2 4 5];
% ax.YLim = [.75 0.95];
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', High Switch vs. No Switch P= ' num2str(round(p12,3)) ', Low Switch vs. No Switch P= ' num2str(round(p13,3))  ]);
%
%
%
%
%
%
%
% figure(1831), axis square, hold on;
% errorbar([1 2],[ nanmean((HV1+HA1)./2-n1h) nanmean((HV2+HA2)./2-n1l) ],...
%     [nanstd((HV1+HA1)./2-n1h) nanstd((HV2+HA2)./2-n1l) ]./sqrt(length(Pp)),'r','linewidth',3), hold on
%
% errorbar([4 5 ],[  nanmean((LV2+LA2)./2-n1l) nanmean((LV1+LA1)./2-n1h)],...
%     [ nanstd((LV2+LA2)./2-n1l)  nanstd((LV1+LA1)./2-n1h)]./sqrt(length(Pp)),'b','linewidth',3)
%
% plot([0 6],[0 0],':k')
%
% ax = gca;
% ax.XTickLabel = {'No-Switch-H','Switch-H','No-Switch-L','Switch-L'};
% ax.XTick = [ 1 2 4 5];
% %ax.YLim = [-0.02 0.02];
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', High Switch vs. No Switch P= ' num2str(round(p120,3)) ', Low Switch vs. No Switch P= ' num2str(round(p130,3))  ]);
%
% figure(1732), axis square, hold on;
% errorbar([1 2],[ nanmean((HV1)-n1h) nanmean((HV2)-n1l) ],...
%     [nanstd((HV1)-n1h) nanstd((HV2)-n1l) ]./sqrt(length(Pp)),'r','linewidth',3), hold on
%
% errorbar([4 5 ],[  nanmean((LV2)-n1l) nanmean((LV1)-n1h)],...
%     [ nanstd((LV2)-n1l)  nanstd((LV1)-n1h)]./sqrt(length(Pp)),'b','linewidth',3)
%
% errorbar([1.1 2.1],[ nanmean((HA1)-n1h) nanmean((HA2)-n1l) ],...
%     [nanstd((HV1)-n1h) nanstd((HV2)-n1l) ]./sqrt(length(Pp)),':r','linewidth',3), hold on
%
% errorbar([4.1 5.1 ],[  nanmean((LA2)-n1l) nanmean((LA1)-n1h)],...
%     [ nanstd((LA2)-n1l)  nanstd((LA1)-n1h)]./sqrt(length(Pp)),':b','linewidth',3)
%
% plot([0 6],[0 0],':k')
%
% ax = gca;
% ax.XTickLabel = {'No-Switch-H','Switch-H','No-Switch-L','Switch-L'};
% ax.XTick = [ 1 2 4 5];
% %ax.YLim = [-0.02 0.02];
% title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
%     ', N= ' num2str(length(Pp)) ', HS_HNS_V P= ' num2str(round(p1200,3)) ', HS_HNS_A P= ' num2str(round(p1201,3))...
%     ', LS_LNS_V P= ' num2str(round(p1300,3)) ', LS_LNS_A P= ' num2str(round(p1301,3))]);


figure(1833), axis square, hold on;
errorbar([4 5],[ nanmean((HV1+HA1)./2-(HV_pre+HA_pre)./2) nanmean((HV2+HA2)./2-(HV_pre+HA_pre)./2) ],...
    [nanstd((HV1+HA1)./2-(HV_pre+HA_pre)./2) nanstd((HV2+HA2)./2-(HV_pre+HA_pre)./2) ]./sqrt(length(Pp)),'r','linewidth',3), hold on
errorbar([1.1 2.1],[ nanmean(HV1-HV_pre) nanmean(HV2-HV_pre) ],...
    [nanstd(HV1-HV_pre) nanstd(HV2-HV_pre) ]./sqrt(length(Pp)),'c','linewidth',3), hold on
errorbar([1.2 2.2],[ nanmean(HA1-HA_pre) nanmean(HA2-HA_pre) ],...
    [nanstd(HA1-HA_pre) nanstd(HA2-HA_pre)]./sqrt(length(Pp)),'m','linewidth',3), hold on

% errorbar([4 5],[  nanmean((LV1+LA1)./2-(LV_pre+LA_pre)./2) nanmean((LV2+LA2)./2-(LV_pre+LA_pre)./2)],...
%     [ nanstd((LV1+LA1)./2-(LV_pre+LA_pre)./2)  nanstd((LV2+LA2)./2-(LV_pre+LA_pre)./2)]./sqrt(length(Pp)),'b','linewidth',3)
% errorbar([4.1 5.1],[ nanmean(LV1-LV_pre) nanmean(LV2-LV_pre) ],...
%     [nanstd(LV1-LV_pre) nanstd(LV2-LV_pre) ]./sqrt(length(Pp)),'c','linewidth',3), hold on
% errorbar([4.2 5.2],[ nanmean(LA1-LA_pre) nanmean(LA2-LA_pre) ],...
%     [nanstd(LA1-LA_pre) nanstd(LA2-LA_pre) ]./sqrt(length(Pp)),'m','linewidth',3), hold on
%
% errorbar([6 7],[  nanmean(N1-N_pre) nanmean(N2-N_pre)],...
%     [ nanstd(N1-N_pre)  nanstd(N2-N_pre)]./sqrt(length(Pp)),'k','linewidth',3)

errorbar([4 5],[  nanmean((LV1+LA1)./2-(LV_pre+LA_pre)./2) nanmean((LV2+LA2)./2-(LV_pre+LA_pre)./2)],...
    [ nanstd((LV1+LA1)./2-(LV_pre+LA_pre)./2)  nanstd((LV2+LA2)./2-(LV_pre+LA_pre)./2)]./sqrt(length(Pp)),'b','linewidth',3)
errorbar([1.1 2.1],[ nanmean(LV1-LV_pre) nanmean(LV2-LV_pre) ],...
    [nanstd(LV1-LV_pre) nanstd(LV2-LV_pre) ]./sqrt(length(Pp)),'c','linewidth',3), hold on
errorbar([1.2 2.2],[ nanmean(LA1-LA_pre) nanmean(LA2-LA_pre) ],...
    [nanstd(LA1-LA_pre) nanstd(LA2-LA_pre) ]./sqrt(length(Pp)),'m','linewidth',3), hold on

errorbar([6 7],[  nanmean(N1-N_pre) nanmean(N2-N_pre)],...
    [ nanstd(N1-N_pre)  nanstd(N2-N_pre)]./sqrt(length(Pp)),'k','linewidth',3)


ax = gca;
ax.XTickLabel = {'No-Switch-H','Switch-H','No-Switch-L','Switch-L','Switch-previousH','Switch-previousL'};
ax.XTickLabel = {'after-H','After-L','After-H','After-L','After-H','After-L'};

ax.XTick = [ 1 2 4 5 6 7];
ax.YLim = [-.08 0.15];
title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
    ', N= ' num2str(length(Pp)) ', High Switch vs. No Switch P= ' num2str(round(p12,3)) ', Low Switch vs. No Switch P= ' num2str(round(p13,3))  ]);
plot([0 6],[0 0],':k')
[h p12000]= ttest((HV1+HA1)./2-(HV_pre+HA_pre)./2)  %%%%% this is switch versus no switch
[h p12001]= ttest((HV2+HA2)./2-(HV_pre+HA_pre)./2)  %%%%% this is switch versus no switch

[h p12002]= ttest((LV1+LA1)./2-(LV_pre+LA_pre)./2)  %%%%% this is switch versus no switch
[h p12003]= ttest((LV2+LA2)./2-(LV_pre+LA_pre)./2)  %%%%% this is switch versus no switch
title([ experiment_name ', Blocks:' prepostVec{prepost} ' ' groupVec{group}...
    ', N= ' num2str(length(Pp)) ', Learning_H-Noswith P= ' num2str(round(p12000,3)) ', Learning_H-swith P= ' num2str(round(p12001,3))...
    ', Learning_L-Noswith P= ' num2str(round(p12002,3)) ', Learning_L-swith  P= ' num2str(round(p12003,3))]);
[h p12004]= ttest(((HV1+HA1)./2-(HV_pre+HA_pre)./2)-((LV2+LA2)./2-(LV_pre+LA_pre)./2 ))
[h p12005]= ttest(((HA1+HA1)./2-(HA_pre+HA_pre)./2)-((LA1+LA1)./2-(LA_pre+LA_pre)./2 ))

%%
%%%%%%%%%%%%%%%%%%%%% to look at all trials and account for the effect of
%%%%%%%%%%%%%%%%%%%%% the previous trial see figure 11833
pre=load('pre_beh5_jessica_RT.mat','s','Pp');
HV_pre= pre.s.subBCvH(Pp)';
LV_pre= pre.s.subBCvL(Pp)';
HA_pre= pre.s.subSPvH(Pp)';
LA_pre= pre.s.subSPvL(Pp)';
N_pre = pre.s.subNeut(Pp)';

% HV_pre= zeros(length(Pp),1);
% LV_pre= zeros(length(Pp),1);
% HA_pre= zeros(length(Pp),1);
% LA_pre= zeros(length(Pp),1);
% N_pre= zeros(length(Pp),1);

maxn = max(max(cellfun(@length,s.vecAccu)));
HV = nan(length(Pp),maxn);
LV = nan(length(Pp),maxn);
HA = nan(length(Pp),maxn);
LA = nan(length(Pp),maxn);
N = nan(length(Pp),maxn);
subj = nan(maxn,length(Pp),4);
rewconds={'H';'L';'H';'L'}';
modalityconds={'Visual';'Visual';'Auditory';'Auditory'}';
binconds={'1st';'2ns'}';
prevtrial = {'high';'low';'neutr'};
prevtrialacu = {'err';'corr'};
%prevtrialacu2 = {'high_ncorrect';'high_correct';'low_ncorrect';'low_correct';'neut'};
prevtrialacu2 = {'high-err';'high-cor';'low-err';'low-cor';'neut'};
prevtrialmod = {'Visual';'Auditory';'Neutr'};
rewfactor = cell(maxn,length(Pp),length(rewconds));
modalityfactor = cell(maxn,length(Pp),length(rewconds));
binfactor = cell(maxn,length(Pp),length(binconds));
accufac = nan(maxn,length(Pp),length(rewconds));
prevtrialfac = cell(maxn,length(Pp),length(binconds));
prevtrialacufac = cell(maxn,length(Pp),length(binconds));
prevtrialacufac2 = cell(maxn,length(Pp),length(binconds));
prevtrialmodfac2 = cell(maxn,length(Pp),length(binconds));

for p=1:length(Pp)
    HV(p,1:length(s.vecAccu{Pp(p),3}))= s.vecRT{Pp(p),3}-HV_pre(p);
    LV(p,1:length(s.vecAccu{Pp(p),2}))= s.vecRT{Pp(p),2}-LV_pre(p);
    HA(p,1:length(s.vecAccu{Pp(p),5}))= s.vecRT{Pp(p),5}-HA_pre(p);
    LA(p,1:length(s.vecAccu{Pp(p),4}))= s.vecRT{Pp(p),4}-LA_pre(p);
    N(p,1:length(s.vecAccu{Pp(p),1})) = s.vecRT{Pp(p),1}-N_pre(p);
    subj(1:maxn,p,1:4)= p;
    for j =1:length(rewconds)
        rewfactor(1:maxn,p,j) = [ repmat(rewconds(j),maxn,1)]; %%%%
        modalityfactor(1:maxn,p,j) = [ repmat(modalityconds(j),maxn,1)]; %%%%
        binfactor(1:maxn,p,j) = [ repmat(binconds(1),floor(maxn/2),1) ; repmat(binconds(2),ceil(maxn/2),1)]; %%%%
        accufac(1:maxn,p,j) = [ repmat(s.subAccu(Pp(p)),maxn,1)];
        %%%%%% this way for high reward trials we use labels according
        %%%%%% to wehtehr reward assignments stayed the same or not. This is as opposed to look at trials before high or low
        if j==1 | j==3
            repvec = [2 1 0];
        elseif j==2 | j==4
            repvec = [2 1 0];
        end
        
        for n =1:maxn
            if n== 1
                prevtrialfac(n,p,j) = prevtrial(3); %%%%
                prevtrialacufac(n,p,j) = prevtrialacu(1); %%%%
                prevtrialacufac2(n,p,j) = prevtrialacu2(5); %%%%
                prevtrialmodfac2 (n,p,j) = prevtrialmod(3); %%%%
            elseif n>1 & n<= length(s.vecAccu{Pp(p),j})
                %%%%%%%%%%%%%%%%%%%%% consider only reward  of previous trial
                if s.vecCondRew0{Pp(p),j}(n)  == repvec(1)
                    prevtrialfac (n,p,j) = prevtrial(1); %%%%
                elseif s.vecCondRew0{Pp(p),j}(n)  == repvec(2)
                    prevtrialfac (n,p,j) = prevtrial(2); %%%%
                elseif s.vecCondRew0{Pp(p),j}(n)  == 0
                    prevtrialfac (n,p,j) = prevtrial(3); %%%%
                end
                %%%%%%%%%%%%%%%%%%%%% consider only accuracy of previous trial
                if s.vecAccu0{Pp(p),j}(n)  == 1
                    prevtrialacufac(n,p,j) = prevtrialacu(2); %%%%
                elseif s.vecAccu0{Pp(p),j}(n)  == 0
                    prevtrialacufac(n,p,j) = prevtrialacu(1); %%%%
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%this is the main thing: consider the conjunction of  accuracy and reward of previous trial
                if s.vecCondRew0{Pp(p),j}(n)  == repvec(1) && s.vecAccu0{Pp(p),j}(n)  == 1
                    prevtrialacufac2 (n,p,j) = prevtrialacu2(2); %%%%
                elseif s.vecCondRew0{Pp(p),j}(n)  == repvec(1) && s.vecAccu0{Pp(p),j}(n)  == 0
                    prevtrialacufac2 (n,p,j) = prevtrialacu2(1); %%%%
                elseif s.vecCondRew0{Pp(p),j}(n)  == repvec(2) && s.vecAccu0{Pp(p),j}(n)  == 1
                    prevtrialacufac2 (n,p,j) = prevtrialacu2(4); %%%%
                elseif s.vecCondRew0{Pp(p),j}(n)  == repvec(2) && s.vecAccu0{Pp(p),j}(n)  == 0
                    prevtrialacufac2 (n,p,j) = prevtrialacu2(3); %%%%
                elseif s.vecCondRew0{Pp(p),j}(n)  == repvec(3)
                    prevtrialacufac2 (n,p,j) = prevtrialacu2(5); %%%%
                end
                
                
                
                if s.vecCondMod0{Pp(p),j}(n)  == repvec(1)
                    prevtrialmodfac2 (n,p,j) = prevtrialmod(2); %%%%
                elseif s.vecCondMod0{Pp(p),j}(n)  == repvec(2)
                    prevtrialmodfac2 (n,p,j) = prevtrialmod(1); %%%%
                elseif s.vecCondMod0{Pp(p),j}(n)  == repvec(3)
                    prevtrialmodfac2 (n,p,j) = prevtrialmod(3); %%%%
                end
                %%%%%%%%%%%%%%%%%%%%%this is the main thing
            elseif n> length(s.vecCondRew0{Pp(p),j})
                prevtrialfac (n,p,j) = prevtrial(3); %%%%
                prevtrialacufac(n,p,j) = prevtrialacu(1); %%%%
                prevtrialacufac2(n,p,j) = prevtrialacu2(5); %%%%
                prevtrialmodfac2 (n,p,j) = prevtrialmod(3); %%%%
            end
        end
        
    end
end
anovdata = cat(3,HV',LV', HA', LA');

%%%%%now reshape them: column vector

anovdata = anovdata(:);
subj = subj(:); %%%%;
rewfactor = rewfactor(:); %%%%
modalityfactor = modalityfactor(:); %%%%
binfactor = binfactor(:); %%%%
prevtrialfac = prevtrialfac(:); %%%%
prevtrialacufac = prevtrialacufac(:);
prevtrialacufac2 = prevtrialacufac2(:);
accufac =  accufac(:);
prevtrialmodfac2 = prevtrialmodfac2(:);
%%%%%%%%%%%%%this one should be worked out
% [h p1]= anovan(anovdata,{subj rewfactor  modalityfactor },'model',[1 0 0  ;0 1 0  ;0 0 1   ;0 1 1   ],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'    }); %%%%% repeated measure ANOVA

% [h p2 stats1]= anovan(anovdata,{subj rewfactor  modalityfactor prevtrialfac },'model',[1 0 0 0 ;0 1 0 0 ;0 0 1 0 ;0 0 0 1 ;0 1 1 0 ;0 1 0 1 ;...
%     0 0 1 1 ;0 1 1 1 ],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'previous'  }); %%%%% repeated measure ANOVA

% [h p3]= anovan(anovdata,{subj rewfactor  modalityfactor prevtrialfac prevtrialacufac },'model',[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;...
%     0 1 1 0 0;0 1 0 1 0;0 1 0 0 1 ;...
%     0 0 1 1 0 ;0 0 1 0 1; 0 1 1 1 1 ],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'previous' 'previousacu' }); %%%%% repeated measure ANOVA
%
% [h p4]= anovan(anovdata,{subj rewfactor  modalityfactor  prevtrialacufac },'model',[1 0 0 0;0 1 0  0;0 0 1  0;0 0 0  1;...
%     0 1 1  0 ; 0 1 0  1 ;...
%     0 0 1 1  ; 0 1 1  1 ],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'previousacu' }); %%%%% repeated measure ANOVA

% [h p22 stats4]= anovan(anovdata,{subj rewfactor  modalityfactor prevtrialfac prevtrialacufac },'model',[1 0 0 0 0 ;0 1 0 0 0;0 0 1 0 0 ;0 0 0 1 0;0 0 0 0 1;...
%     0 1 1 0 0 ;0 1 0 1 0; 0 1 0 0 1;...
%     0 0 1 1 0; 0 0 1 0 1 ; 0 1 0 1 1; 0 0 1 1 1;0 1 1 1 1 ],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'moalityfactor' 'previus'  'PreviousAcu'  }); %%%%% repeated measure ANOVA


% [h p3]= anovan(anovdata,{subj rewfactor  modalityfactor prevtrialfac binfactor},'model',[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0 ;0 0 0 0 1; ...
%     0 1 1 0 0;0 1 0 1 0;0 1 0 0 1;...
%     0 0 1 1 0;0 0 1 0 1;...
%     0 1 1 1 1],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'previous'  'bins'}); %%%%% repeated measure ANOVA
%
% [h p4]= anovan(anovdata,{subj rewfactor  modalityfactor binfactor },'model',[1 0 0 0 ;0 1 0 0 ;0 0 1 0 ;0 0 0 1 ;0 1 1 0 ;0 1 0 1 ;...
%     0 0 1 1 ;0 1 1 1 ],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor' 'bins'  }); %%%%% repeated measure ANOVA

tbl2 = table(subj,prevtrialfac,prevtrialacufac,prevtrialacufac2,prevtrialmodfac2, binfactor, modalityfactor,rewfactor,accufac,anovdata,'VariableNames',{'subjs' 'previous' 'PreviousAcu' 'PreviousAcuRew' 'PreviousMod' 'bins' 'modalityfactor','rewardfactor','genaccu','accuracy' });
%%%%%%%%%%%%%%% to set a certain order for dummyfactors codes: you can
%%%%%%%%%%%%%%% comment this
tbl2.PreviousAcuRew = categorical(tbl2.PreviousAcuRew,{'high-err';'high-cor';'low-err';'low-cor';'neut'});
tbl2.modalityfactor = categorical(tbl2.modalityfactor,{'Visual';'Auditory'}); %%%% visual first and then auditory see: https://de.mathworks.com/matlabcentral/answers/146623-reference-dummy-coding-with-matlab-fitlme
tbl2.rewardfactor = categorical(tbl2.rewardfactor,{'H';'L'});
tbl2.PreviousMod = categorical(tbl2.PreviousMod,{'Visual';'Auditory';'Neutr'}); %%%% visual first and then auditory see: https://de.mathworks.com/matlabcentral/answers/146623-reference-dummy-coding-with-matlab-fitlme
tbl2.previous = categorical(tbl2.previous,{'high';'low';'neutr'}); %%%% visual first and then auditory see: https://de.mathworks.com/matlabcentral/answers/146623-reference-dummy-coding-with-matlab-fitlme

% lme0 = fitlme(tbl2,'accuracy~genaccu+(rewardfactor|subjects)+(modalityfactor|subjects)+(1|subjects)')
% lme1 = fitlme(tbl2,'accuracy~1+(rewardfactor|subjs)+(modalityfactor|subjs)')%%%%%%%%%%% only random effect has nans because data has nans
% lme2 = fitlme(tbl2,'accuracy~1+(rewardfactor|subjs)+(modalityfactor|subjs)+ (previous|subjs)')%%%%%%%%%%% fuxed and random effects
%lme3 = fitlme(tbl2,'accuracy~1+(rewardfactor+modalityfactor+PreviousAcuRew)+(1|subjs)')%%%%%%%%%%% fuxed and random effects
% lme4 = fitlme(tbl2,'accuracy~1+(rewardfactor|subjs)+(modalityfactor|subjs)+ (PreviousAcuRew|subjs)')%%%%%%%%%%% fuxed and random effects
%% Comparison of the results of ANOVAN and linear mixed model (which should be identical)
%%%%%%%%%%%%%%%without interaction term: they re identical
[h p22 stats4]= anovan(anovdata,{subj rewfactor  modalityfactor prevtrialacufac2 },'model',[1 0 0 0 ;0 1 0 0 ;0 0 1 0 ;0 0 0 1 ],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'PreviousAcuRew'  }) %%%%% repeated measure ANOVA
lme3 = fitlme(tbl2,'accuracy~1+(rewardfactor+modalityfactor+PreviousAcuRew)+(1|subjs)','DummyVarCoding','effects')%%%%%%%%%%% fuxed and random effects

%%%%%%%%%%%%%%%with interaction term
%%%%Comparison of the results of ANOVAN and linear mixed model (which should be identical)
%%%%%%%%%%%%%%%with interaction term
%%%%%%%%%%%%%%%% the model assumes main and interaction effect of reward
%%%%%%%%%%%%%%%% and modality and enters previoustrial factors as a
%%%%%%%%%%%%%%%% covariate too
[h p220 stats40]= anovan(anovdata,{subj rewfactor  modalityfactor prevtrialacufac2 },'model',[1 0 0 0 ;0 1 0 0 ;0 0 1 0 ;0 0 0 1 ;0 1 1 0  ],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'PreviousAcuRew'  }) %%%%% repeated measure ANOVA
lme4 = fitlme(tbl2,'accuracy~(rewardfactor*modalityfactor+PreviousAcuRew)+ (1|subjs)','DummyVarCoding','effects')%%%%%%%%%%% fuxed and random effects: check DummyVArCoding
%%%%%%%%%%%%%%%% main reference: https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
%%%%%%%%%%%%https://revistas.usb.edu.co/index.php/IJPR/article/view/807
%%%%%%%%%%http://www.bodowinter.com/tutorial/bw_LME_tutorial2.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% very good review https://psych.wisc.edu/Brauer/BrauerLab/wp-content/uploads/2014/04/Brauer_and_Curtin_LMEMs-2017-Psych_Methods.pdf
%%%%%%%%%%%%%%%%%https://www.frontiersin.org/articles/10.3389/fpsyg.2015.00002/full
%%%%%%%%%%%%%https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5366825/#
%%%%%%%%%%%see https://stats.stackexchange.com/questions/282646/two-way-repeated-measures-linear-mixed-model
%%%%%%%%%%%%%%%%%%to visualize this
%%%%%%%%%%%%%%%%%%%%%https://stats.stackexchange.com/questions/387699/which-random-effects-to-include-in-a-mixed-effects-model
figure(1900),bar([1 2 ],[nanmean(nanmean(HV,2)) nanmean(nanmean(LV,2)) ],'c'), hold on
plot([1 2 ],[(nanmean(HV,2)) (nanmean(LV,2)) ],'oc'), hold on
bar([ 4 5],[ nanmean(nanmean(HA,2)) nanmean(nanmean(LA,2))],'m'), hold on
plot([4 5 ],[(nanmean(HA,2)) (nanmean(LA,2)) ],'om'), hold on
ylabel('Performance difference from base line')
%%%%%%%%%%%%%%%% to plot fitlme estimates
figure(2000),bar(1:7,double(lme4.Coefficients(2:end,2)),'y'), hold on
plot([1:7; 1:7],[double(lme4.Coefficients(2:end,end-1)), double(lme4.Coefficients(2:end,end))]','k'), hold on
set (gca, 'xtick', 1:7, 'xticklabel', lme4.CoefficientNames(2:end))

% % %%%%%%%%%%%if we want to include other fixed effects (e.g. between subjects factor such as overall accuracy of each subject) for which the subjs
% % %%%%%%%%%%%should not be considered as random slope we do
% %  %%%%%%%%%%this is based on the example from Matlab help on lme
% % %  Fit a linear mixed-effects model for miles per gallon in the city, with fixed effects for horsepower, and uncorrelated random effect for intercept and horsepower grouped by the engine type.
% % %
% % % lme = fitlme(tbl,'CityMPG~Horsepower+(1|EngineType)+(Horsepower-1|EngineType)');
% %
% %
% %  lme5 = fitlme(tbl2,'accuracy~(rewardfactor*modalityfactor+PreviousAcuRew+genaccu)+ (1|subjs)+(genaccu-1|subjs)','DummyVarCoding','effects')%%%%%%%%%%% fuxed and random effects: check DummyVArCoding

% lme6 = fitlme(tbl2,'accuracy~(rewardfactor*modalityfactor+PreviousAcuRew)+ (rewardfactor*modalityfactor+PreviousAcuRew|subjs)','DummyVarCoding','effects')%%%%%%%%%%% fixed and random effects: check DummyVArCoding
% lme6 = fitlme(tbl2,'accuracy~(1+rewardfactor*modalityfactor+PreviousAcuRew)+ (1+rewardfactor*modalityfactor+PreviousAcuRew|subjs)','DummyVarCoding','effects')%%%%%%%%%%% fixed and random effects: check DummyVArCoding
lme5 = fitlme(tbl2,'accuracy~(rewardfactor*PreviousAcuRew+modalityfactor)+ (1|subjs)','DummyVarCoding','effects')%%%%%%%%%%%
lme6 = fitlme(tbl2,'accuracy~(rewardfactor*modalityfactor*PreviousAcuRew)+ (1|subjs)','DummyVarCoding','effects')%%%%%%%%%%% fixed and random effects: check DummyVArCoding:https://hal.archives-ouvertes.fr/hal-01668131/document
lme7 = fitlme(tbl2,'accuracy~(rewardfactor+PreviousAcuRew+modalityfactor)+ (rewardfactor+PreviousAcuRew+modalityfactor|subjs)','DummyVarCoding','effects')%%%%%%%%%%%
lme9 = fitlme(tbl2,'accuracy~(rewardfactor*modalityfactor)+ (rewardfactor*PreviousAcuRew*modalityfactor|subjs)','DummyVarCoding','effects')%%%%%%%%%%%

lme8 = fitlme(tbl2,'accuracy~(rewardfactor*PreviousAcuRew*modalityfactor)+ (rewardfactor*PreviousAcuRew*modalityfactor|subjs)','DummyVarCoding','effects')%%%%%%%%%%% https://www.frontiersin.org/articles/10.3389/fnhum.2017.00491/full

bla=compare(lme6,lme8)%%%%%%%%%%%%%%% model 6 and 8 are the best models and model 6 could be selected with only random intercept see https://www.nature.com/articles/s41598-018-31526-y?WT.feed_name=subjects_perception
%%%%%%%%%%%%%%%Random slopes for within-subject factors (trial type) were dropped based on likelihood ratio tests on nested models

lme16 = fitlme(tbl2,'accuracy~(rewardfactor*modalityfactor*previous*PreviousAcu)+ (1|subjs)','DummyVarCoding','effects')%%%%%%%%%%% fixed and random effects: check DummyVArCoding
% lme160 = fitlme(tbl2,'accuracy~(rewardfactor*modalityfactor*previous*PreviousAcu)+ (rewardfactor*modalityfactor*previous*PreviousAcu|subjs)','DummyVarCoding','effects')%%%%%%%%%%% fixed and random effects: check DummyVArCoding
% lme1600 = fitlme(tbl2,'accuracy~(rewardfactor*modalityfactor*previous)+ (rewardfactor*modalityfactor*previous*PreviousAcu|subjs)','DummyVarCoding','effects')%%%%%%%%%%% fixed and random effects: check DummyVArCoding

%%
%%%%%%%%%%%%%%%% to plot fitlme estimates
nb =size(lme160.Coefficients,1)-1;
mylm =lme160.Coefficients;
mylm_names = lme160.CoefficientNames;
figure(2000),bar(1:nb,double( mylm(2:end,2)),'y'), hold on
plot([1:nb; 1:nb],[double( mylm(2:end,end-1)), double( mylm(2:end,end))]','k'), hold on
set (gca, 'xtick', 1:nb, 'xticklabel', mylm_names(2:end))

%%
% %%%%%%%%%%%%%%%%%%%%% to look at all trials and account for the effect of
% %%%%%%%%%%%%%%%%%%%%% the previous trial see figure 11833: Not sure
% pre=load('pre_beh5_jessica.mat','s','Pp');
% HV_pre= pre.s.subBCvH(Pp)';
% LV_pre= pre.s.subBCvL(Pp)';
% HA_pre= pre.s.subSPvH(Pp)';
% LA_pre= pre.s.subSPvL(Pp)';
% N_pre= pre.s.subNeut(Pp)';
%
% % HV_pre= zeros(length(Pp),1);
% % LV_pre= zeros(length(Pp),1);
% % HA_pre= zeros(length(Pp),1);
% % LA_pre= zeros(length(Pp),1);
% % N_pre= zeros(length(Pp),1);
%
% maxn = max(max(cellfun(@length,s.vecAccu)));
% HV = nan(length(Pp),maxn);
% LV = nan(length(Pp),maxn);
% HA = nan(length(Pp),maxn);
% LA = nan(length(Pp),maxn);
% N = nan(length(Pp),maxn);
% subj = nan(maxn,length(Pp),4);
% rewconds={'H';'L';'H';'L'}';
% modalityconds={'Visual';'Visual';'Auditory';'Auditory'}';
% binconds={'1st';'2ns'}';
% prevtrial = {'high';'low';'neutr'};
% rewfactor = cell(maxn,length(Pp),length(rewconds));
% modalityfactor = cell(maxn,length(Pp),length(rewconds));
% binfactor = cell(maxn,length(Pp),length(binconds));
% prevtrialfac = cell(maxn,length(Pp),length(binconds));
% accufac = nan(maxn,length(Pp),length(rewconds));
% boxi = cell(maxn,length(Pp),length(binconds));
% pitchi = cell(maxn,length(Pp),length(binconds));
%
% %s.subRewS(Pp,1) s.subRewV(Pp,1)
% for p=1:length(Pp)
%     HV(p,1:length(s.vecAccu{Pp(p),3}))= s.vecAccu{Pp(p),3}-HV_pre(p);
%     LV(p,1:length(s.vecAccu{Pp(p),2}))= s.vecAccu{Pp(p),2}-LV_pre(p);
%     HA(p,1:length(s.vecAccu{Pp(p),5}))= s.vecAccu{Pp(p),5}-HA_pre(p);
%     LA(p,1:length(s.vecAccu{Pp(p),4}))=s.vecAccu{Pp(p),4}-LA_pre(p);
%     N(p,1:length(s.vecAccu{Pp(p),1}))= s.vecAccu{Pp(p),1}-N_pre(p);
%     subj(1:maxn,p,1:4)= p;
%     for j =1:length(rewconds)
%         rewfactor(1:maxn,p,j) = [ repmat(rewconds(j),maxn,1)]; %%%%
%         modalityfactor(1:maxn,p,j) = [ repmat(modalityconds(j),maxn,1)]; %%%%
%         binfactor(1:maxn,p,j) = [ repmat(binconds(1),floor(maxn/2),1) ; repmat(binconds(2),ceil(maxn/2),1)]; %%%%
%         accufac(1:maxn,p,j) = [ repmat(s.subAccu(Pp(p)),maxn,1)]; %%%%
%         %%%%%% this way for high reward trials we use labels according
%         %%%%%% to wehtehr reward assignments stayed the same or not. This is as opposed to look at trials before high or low
%         if j==1 | j==3
%             repvec = [2 1 0];
%         elseif j==2 | j==4
%             repvec = [1 2 0];
%         end
%
%         for n =1:maxn
%             if n== 1
%                 prevtrialfac(n,p,j) = prevtrial(3); %%%%
%             elseif n<= length(s.vecAccu{Pp(p),j})
%                 if s.vecCondRew0{Pp(p),j}(n)  == repvec(1)
%                     prevtrialfac (n,p,j) = prevtrial(1); %%%%
%                 elseif s.vecCondRew0{Pp(p),j}(n)  == repvec(2)
%                     prevtrialfac (n,p,j) = prevtrial(2); %%%%
%                 elseif s.vecCondRew0{Pp(p),j}(n)  == 0
%                     prevtrialfac (n,p,j) = prevtrial(3); %%%%
%                 end
%             elseif n> length(s.vecCondRew0{Pp(p),j})
%                 prevtrialfac (n,p,j) = prevtrial(3); %%%%
%             end
%         end
%
%     end
% end
% anovdata = cat(3,HV',LV', HA', LA');
%
% %%%%%now reshape them: column vector
%
% anovdata = anovdata(:);
% subj = subj(:); %%%%;
% rewfactor = rewfactor(:); %%%%
% modalityfactor = modalityfactor(:); %%%%
% binfactor = binfactor(:); %%%%
% prevtrialfac = prevtrialfac(:); %%%%
%  accufac =  accufac(:);
% % %%%%%%%%%%%%%this one should be worked out
% % [h p]= anovan(anovdata,{subj rewfactor  modalityfactor prevtrialfac },'model',[1 0 0 0 ;0 1 0 0 ;0 0 1 0 ;0 0 0 1 ;0 1 1 0 ;0 1 0 1 ;...
% %     0 0 1 1 ;0 1 1 1 ],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'previous'  }); %%%%% repeated measure ANOVA
% %
% % [h p]= anovan(anovdata,{subj rewfactor  modalityfactor prevtrialfac binfactor},'model',[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0 ;0 0 0 0 1; ...
% %     0 1 1 0 0;0 1 0 1 0;0 1 0 0 1;...
% %     0 0 1 1 0;0 0 1 0 1;...
% %     0 1 1 1 1],'random',1,'varnames',{ 'Subject' 'rewardfactor'   'modalityfactor'  'previous'  'bins'}); %%%%% repeated measure ANOVA
% tbl1 = table(prevtrialfac,modalityfactor,rewfactor,anovdata,'VariableNames',{ 'previous', 'modalityfactor','rewardfactor','accuracy' });
% tbl2 = table(subj,prevtrialfac,modalityfactor,rewfactor,accufac,anovdata,'VariableNames',{'subjects' 'previous' 'modalityfactor','rewardfactor','genaccu','accuracy' });
%
% %lm = fitlm(tbl1,'accuracy~1+rewardfactor*modalityfactor+modalityfactor*previous+rewardfactor*previous+rewardfactor*modalityfactor*previous')
% lm = fitlm(tbl1,'accuracy~1+rewardfactor*modalityfactor*previous')
%
% tbl = anova(lm,'summary')
% %lm = fitlm(tbl2,'accuracy~1+rewardfactor+modalityfactor')
% lm = fitlm(tbl2,'accuracy~1+rewardfactor*modalityfactor')
%
%
% lme0 = fitlme(tbl2,'accuracy~genaccu+(rewardfactor|subjects)+(modalityfactor|subjects)+(1|subjects)')
% lme1 = fitlme(tbl2,'accuracy~1+(1+rewardfactor|subjects)+(1+modalityfactor|subjects)')%%%%%%%%%%% only random effect has nans because data has nans
% lme2 = fitlme(tbl2,'accuracy~rewardfactor*modalityfactor+(1+rewardfactor|subjects)+(1+modalityfactor|subjects)+(1|subjects)')%%%%%%%%%%% fuxed and random effects
%
% [~,~,STATS] = randomEffects(lme1); % Compute the random-effects statistics (STATS)
% % STATS.Level = nominal(STATS.Level);
% % K = zeros(length(STATS),1);
% % K(STATS.Level(tbl2.rewardfactor(:)=='H')) = 1;
% % pVal = coefTest(lme,[0 0 0 0 0 0 0 0 0],0,'REContrast',K')
% %%
% pre=load('pre_beh5_jessica.mat','s','Pp');
% post1H=load('post_beh5_jessica_afterHigh_nminus1_error.mat','s','Pp');
% post2H=load('post_beh5_jessica_afterHigh2_error.mat','s','Pp');
%
% post1L=load('post_beh5_jessica_afterLow_nminus1_error.mat','s','Pp');
% post2L=load('post_beh5_jessica_afterLow2_error.mat','s','Pp');
%
%
%
%
% HV1_err= post2H.s.subBCvH(Pp)';
% LV1_err= post2H.s.subBCvL(Pp)';
% HA1_err= post2H.s.subSPvH(Pp)';
% LA1_err= post2H.s.subSPvL(Pp)';
% N1_err= post2H.s.subNeut(Pp)';
%
% HV2_err= post2L.s.subBCvH(Pp)';
% LV2_err= post2L.s.subBCvL(Pp)';
% HA2_err= post2L.s.subSPvH(Pp)';
% LA2_err= post2L.s.subSPvL(Pp)';
% N2_err= post2L.s.subNeut(Pp)';
%
%
% HV_pre= pre.s.subBCvH(Pp)';
% LV_pre= pre.s.subBCvL(Pp)';
% HA_pre= pre.s.subSPvH(Pp)';
% LA_pre= pre.s.subSPvL(Pp)';
% N_pre= pre.s.subNeut(Pp)';
%
%
%
% figure(1720),  hold on;
% bar([1:4 10:13],[ nanmean((HV1+HV2)./2) nanmean((LV1+LV2)./2) nanmean((HA1+HA2)./2) nanmean((LA1+LA2)./2) ...
%    nanmean((HV1_err+HV2_err)./2) nanmean((LV1_err+LV2_err)./2) nanmean((HA1_err+HA2_err)./2) nanmean((LA1_err+LA2_err)./2) ],'w'), axis square, hold on
% % errorbar([1 3 6 8 11 13 16 18 21 23  ],[nanmean(n1h) nanmean(n1l) nanmean(HV1) nanmean(HV2) nanmean(LV2) nanmean(LV1) nanmean(HA1)...
% %      nanmean(HA2) nanmean(LA2)  nanmean(LA1) ],...
% %     [nanstd(n1h) nanstd(n1l) nanstd(HV1) nanstd(HV2) nanstd(LV2) nanstd(LV1) nanstd(HA1)...
% %      nanstd(HA2) nanstd(LA2)  nanstd(LA1) ]./sqrt(length(Pp)),'k','linewidth',3, 'linestyle','none')
% % ax = gca;
% % ax.XTick = [1 3 6 8 11 13 16 18 21 23  ];
% % ax.XTickLabel = {'n-1h','n-1l','VHNS','VHS','VLNS','VLS','AHNS','AHS','ALNS','ALS'};
% % ax.YLim = [.6 1];
% [h p1]= ttest((HV1_err+HV2_err+HA1_err+HA2_err)./4-(LV1_err+LV2_err+LA1_err+LA2_err)./4,(HV1+HV2+HA1+HA2)./4-(LV1+LV2+LA1+LA2)./4)
%
% %%
% % pre=load('pre_beh5_jessica_RT.mat','s','Pp');
% % post1H=load('post_beh5_jessica_afterHigh_nminus1_RT.mat','s','Pp');
% % post2HV=load('post_beh5_jessica_afterHigh2_PostVisual_RT.mat','s','Pp');
% % post2HA=load('post_beh5_jessica_afterHigh2_PostAuditory_RT.mat','s','Pp');
% %
% % post1L=load('post_beh5_jessica_afterLow_nminus1_RT.mat','s','Pp');
% % post2LV=load('post_beh5_jessica_afterLow2_PostVisual_RT.mat','s','Pp');
% % post2LA=load('post_beh5_jessica_afterLow2_PostAuditory_RT.mat','s','Pp');
% %
% %
% % n1h= (post1H.s.subBCvH(Pp)'+post1H.s.subBCvL(Pp)'+post1H.s.subSPvH(Pp)'+ post1H.s.subSPvL(Pp)'+post1H.s.subNeut(Pp)')/5;
% %
% % n1l= (post1L.s.subBCvH(Pp)'+post1L.s.subBCvL(Pp)'+post1L.s.subSPvH(Pp)'+ post1L.s.subSPvL(Pp)'+post1L.s.subNeut(Pp)')/5;
% %
% %
% % HV1_V= post2HV.s.subBCvH(Pp)';
% % LV1_V= post2HV.s.subBCvL(Pp)';
% % HA1_V= post2HV.s.subSPvH(Pp)';
% % LA1_V= post2HV.s.subSPvL(Pp)';
% % N1_V= post2HV.s.subNeut(Pp)';
% %
% %
% % HV1_A= post2HA.s.subBCvH(Pp)';
% % LV1_A= post2HA.s.subBCvL(Pp)';
% % HA1_A= post2HA.s.subSPvH(Pp)';
% % LA1_A= post2HA.s.subSPvL(Pp)';
% % N1_A= post2HA.s.subNeut(Pp)';
% %
% % HV2_V= post2LV.s.subBCvH(Pp)';
% % LV2_V= post2LV.s.subBCvL(Pp)';
% % HA2_V= post2LV.s.subSPvH(Pp)';
% % LA2_V= post2LV.s.subSPvL(Pp)';
% % N2_V= post2LV.s.subNeut(Pp)';
% %
%
% HV2_A= post2LA.s.subBCvH(Pp)';
% LV2_A= post2LA.s.subBCvL(Pp)';
% HA2_A= post2LA.s.subSPvH(Pp)';
% LA2_A= post2LA.s.subSPvL(Pp)';
% N2_A= post2LA.s.subNeut(Pp)';
%
%
% HV_pre= pre.s.subBCvH(Pp)';
% LV_pre= pre.s.subBCvL(Pp)';
% HA_pre= pre.s.subSPvH(Pp)';
% LA_pre= pre.s.subSPvL(Pp)';
% N_pre= pre.s.subNeut(Pp)';
%
% figure(1800), hold on;
% plot([4:7],[nanmean(HV1_V) nanmean(HV2_V) nanmean(HV1_A) nanmean(HV2_A) ],'r','marker','s'), axis square, hold on
% plot([4:7],[nanmean(LV2_V) nanmean(LV1_V) nanmean(LV2_A) nanmean(LV1_A) ],':r','marker','s'), axis square, hold on
% plot([4:7],[nanmean(HA1_A) nanmean(HA2_A) nanmean(HA1_V) nanmean(HA2_V) ],'b','marker','s'), axis square, hold on
% plot([4:7],[nanmean(LA2_A) nanmean(LA1_A) nanmean(LA2_V) nanmean(LA1_V) ],':b','marker','s'), axis square, hold on
%
%
% figure(1801), hold on;
% plot([4:7],[nanmean(HV1_V-LV2_V) nanmean(HV2_V-LV1_V) nanmean(HV1_A-LV2_A) nanmean(HV2_A-LV1_A) ],'r','marker','s'), axis square, hold on
% plot([4:7],[nanmean(HA1_A-LA2_A) nanmean(HA2_A-LA1_A) nanmean(HA1_V-LA2_V) nanmean(HA2_V-LA1_V) ],'b','marker','s'), axis square, hold on
% ax = gca;
%
% ax.XTick = [4:7 ];
% ax.XTickLabel = {'SRSM','DRSM','SRDM','DRDM'};
%
%
% figure(1802), hold on;
% plot([4:7],[nanmean((HV1_V+LV2_V)/2) nanmean((HV2_V+LV1_V)/2) nanmean((HV1_A+LV2_A)/2) nanmean((HV2_A+LV1_A)/2) ],'r','marker','s'), axis square, hold on
% plot([4:7],[nanmean((HA1_A+LA2_A)/2) nanmean((HA2_A+LA1_A)/2) nanmean((HA1_V+LA2_V)/2) nanmean((HA2_V+LA1_V)/2) ],'b','marker','s'), axis square, hold on
%
% ax = gca;
%
% ax.XTick = [4:7 ];
% ax.XTickLabel = {'prevH-SM','prevL-SM','prevH-DM','prevL-DM'};
%
%
% figure(1803), hold on;
% plot([4:7],[nanmean((HV1_V+LV1_V)/2) nanmean((HV2_V+LV2_V)/2) nanmean((HV1_A+LV1_A)/2) nanmean((HV2_A+LV2_A)/2) ],'r','marker','s'), axis square, hold on
% plot([4:7],[nanmean((HA1_A+LA1_A)/2) nanmean((HA2_A+LA2_A)/2) nanmean((HA1_V+LA1_V)/2) nanmean((HA2_V+LA2_V)/2) ],'b','marker','s'), axis square, hold on
%
% ax = gca;
%
% ax.XTick = [4:7 ];
% ax.XTickLabel = {'prevH-SM','prevL-SM','prevH-DM','prevL-DM'};
