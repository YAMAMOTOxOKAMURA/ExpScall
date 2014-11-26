%
%  Play Size disrimination Exp. & GCFB simulation
%  Irino, T.
%  Created:  23 Sep 2014 
%  Modified: 23 Sep 2014 
%  Modified: 24 Sep 2014 
%  Modified: 25 Sep 2014 (IT)  % text manupilation
%  Modified: 27 Sep 2014 (IT)  
%  Modified: 28 Sep 2014 (IT)  % SIMparam.SwStatJudge = 4; for variation 
%  Modified: 29 Sep 2014 (IT)  % StatJudge. ...
%  Modified: 1 Nov 2014 (KY) % SIMparam.SwStatJudge = 6,vowel = /a/
%  Modified: 2 Nov 2014 (KY) % SIMparam.SwStatJudge = 6,vowel = /a-o/
%
%
function [RsltSim] = PlayExpWordSize_Simulation(SndSim,SIMparam);

disp('test');
    if isfield(SIMparam,'Method') == 0, 
       error('Specify SIMparam.Method');
    end;
    if isfield(SIMparam,'fs') == 0, 
       error('Specify SIMparam.fs');
    end;
    disp(sprintf('SIMparam.SwStatJudge   = %g',  SIMparam.SwStatJudge));
       
    fs = SIMparam.fs;
    %% Method
    if SIMparam.Method == 1,
        GCparam.fs    = fs;
        GCparam.NumCh = 100;
        GCparam.OutMidCrct = 'ELC';
        GCparam.Ctrl   = 'dynamic';
        SGparam.Method = 1; % method 1
        GCparam.FRange = [100, 6000];
    else
        error('Not prepared yet. SIMparam.Method');
    end;
    
    SIMparam.ExpSndLeveldB = 70;
    SIMparam.ExdSndDigitalLeveldB = -26;
    [dummy AmpdB] = Eqlz2MeddisHCLevel(1,0);
    SIMparam.AmpdB_MHCL0dB = AmpdB(3);
    SIMparam.LevelCnvtdB = SIMparam.ExpSndLeveldB - SIMparam.AmpdB_MHCL0dB - SIMparam.ExdSndDigitalLeveldB;
    % Digital Level --> MeddisHCLevel 66 dB = 70 - 30 -(-26)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation    %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%    
    NameSpec = [SIMparam.DirSpec  SIMparam.NameSpec ];    
    NameSpecMat = [ NameSpec '.mat'];
    if exist(NameSpecMat) ~= 0,
        
        str = ['load ' NameSpecMat ];
        disp(str);
        eval(str);

    else
      %% Producing spectrogram    
      for nIntvl = 1:2,
        SndSim1 = SndSim(nIntvl,:);
        Snd = 10^(SIMparam.LevelCnvtdB/20)*SndSim1;
        [dummy AmpdB] = Eqlz2MeddisHCLevel(Snd,0);
        disp(['Snd_MeddisHCLevel = ' num2str(AmpdB(3)) ' (dB)']);
        Tsnd = length(Snd)/fs;
        
        tic;
        [cGCout, pGCout, GCparam, GCresp] = GCFBv209(Snd,GCparam);
        %        GCresp.Fr1
        tm = toc;
        disp(['Elapsed time is ' num2str(tm,4) ' (sec) = ' ...
             num2str(tm/Tsnd,4) ' times RealTime.']);
        disp(' ');

        %% Spectrogram
        SGparam.fs = fs;
        [GCSpecGramPwr, SGparam] = CalSmoothSpec(cGCout.^2,SGparam);
        %[GCSpecGramPwr, SGparam] = CalSmoothSpec(pGCout.^2,SGparam);
        [NumCh NumFrame] = size(GCSpecGramPwr);
        GCSpecGram(1:NumCh,1:NumFrame,nIntvl) = sqrt(GCSpecGramPwr); %  rms spectrogram
            
      end;
      str = ['save ' NameSpec ' GCSpecGram SGparam'];
      disp(str);
      eval(str);
       
    end;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%
    % Name extraction    %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    for nIntvl = 1:2 % 2 interval 
      for nRpt = 1:SIMparam.NumRptWord,
        nwd = (nIntvl-1)*SIMparam.NumRptWord + nRpt;
        NameSndPlay = char(SIMparam.NameSndPlayCell(nwd));
        Name_SndPlay(nwd) = cellstr(NameSndPlay(1:8)); % FW03 name
        NumVct = Name2NumFW03(NameSndPlay(1:8));
        [dummy,fs1, dummy3, TextWord] = LoadFW03(NumVct(1),NumVct(2),NumVct(3),NumVct(4),-1);
        TextWord_SndPlay(nwd) = cellstr(TextWord);
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add vowel Estimate (load HTK file) %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        [HTKdata,htktext,SndInf,Jparam] = LoadFW03label(SIMparam,NameSndPlay);
        str = ['Interval ' num2str(nIntvl) '-' num2str(nRpt) ': ' htktext '  Snd_time(sec): ' num2str(SndInf.tsec)];
        disp(str);
        %disp(SndTm);
        %htkorg{1}
        
        Snd_Tsec(nRpt) = SndInf.nFlame*Jparam.Tshift; %%
        Snd_Flame(nRpt) = SndInf.nFlame; 
        vowellabel(1:length(HTKdata),nRpt,nIntvl) = horzcat(HTKdata.type);
        vowelCTsec(1:length(HTKdata),nRpt,nIntvl) = horzcat(HTKdata.center);
        vowelStrtTsec(1:length(HTKdata),nRpt,nIntvl) = horzcat(HTKdata.begin);
        vowelEdTsec(1:length(HTKdata),nRpt,nIntvl) = horzcat(HTKdata.end);
        
        if nRpt == 2,
            vowelCTsec(1:length(HTKdata),nRpt,nIntvl) = vowelCTsec(1:length(HTKdata),nRpt,nIntvl) + Snd_Tsec(1);
            vowelStrtTsec(1:length(HTKdata),nRpt,nIntvl) = vowelStrtTsec(1:length(HTKdata),nRpt,nIntvl) + Snd_Tsec(1);
            vowelEdTsec(1:length(HTKdata),nRpt,nIntvl) = vowelEdTsec(1:length(HTKdata),nRpt,nIntvl) + Snd_Tsec(1);
        end;
        
      end;
      str = ['Interval ' num2str(nIntvl) ' : time is ' num2str(sum(Snd_Tsec))];
      disp(str);
    end;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot Spectrogram  SIMparam.SwPlot = 0 or 1 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    SIMparam.SwPlot = 0;
    
    if SIMparam.SwPlot == 1,
    [NumCh NumFrame NumIntrvl] = size(GCSpecGram);
    subplot(2,1,1);
    MaxValue = max(max(max(GCSpecGram)));
    tsec = (0:NumFrame-1)*SGparam.Tshift;
    image(tsec, 1:NumCh,GCSpecGram(:,:,1)/MaxValue*64);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %tsec = (bigf_1:endf_1)*SGparam.Tshift;
    %image(tsec, 1:NumCh,GCSpecGram(:,bigf_1:endf_1,1)/MaxValue*64);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(gca,'YDir','normal');
    xlabel('Time (sec)'); 
    ylabel('channel');
    title([char(Name_SndPlay(1)) ' "' char(TextWord_SndPlay(1)) '" ---  ' ...
           char(Name_SndPlay(2)) ' "' char(TextWord_SndPlay(2)) '"'], ...
           'Interpreter','none');
       
    subplot(2,1,2);
    image(tsec, 1:NumCh,GCSpecGram(:,:,2)/MaxValue*64);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %tsec = (bigf_2:endf_2)*SGparam.Tshift;
    %image(tsec, 1:NumCh,GCSpecGram(:,bigf_2:endf_2,1)/MaxValue*64);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(gca,'YDir','normal');
    xlabel('Time (sec)');
    ylabel('channel');
    title([char(Name_SndPlay(3)) ' "' char(TextWord_SndPlay(3)) '" ---  ' ...
           char(Name_SndPlay(4)) ' "' char(TextWord_SndPlay(4)) '"'], ...
           'Interpreter','none');       
    str = ['print -depsc ' NameSpec];
    disp(str);
    eval(str);
    drawnow
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Judgement Statistics %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %GCSpecGram
    [NumCh NumFrame NumIntvl] = size(GCSpecGram);
    
    %pause;
    
      %% Stat1: Simple CoG calculation 
    for nIntvl = 1:NumIntvl,
        ExctPtrnAll= rms(GCSpecGram(:,:,nIntvl),2);
        CumExtPtrnAll = cumsum(ExctPtrnAll);
        ChIntrp = 0.01;
        NchRes = 1:ChIntrp:NumCh;
        IntrpCumExctPtrnAll = interp1(1:NumCh,CumExtPtrnAll,NchRes);
        [dummy Ncog1] = min(abs(IntrpCumExctPtrnAll-CumExtPtrnAll(end)/2));
        StatJudge.ChCoG4AllDur(nIntvl) = NchRes(Ncog1);
        
        %%
        
        % power base  28 Sep 14
        CumExtPtrnAll2 = cumsum(ExctPtrnAll.^2);
        IntrpCumExctPtrnAll2 = interp1(1:NumCh,CumExtPtrnAll2,NchRes);
        [dummy Ncog2] = min(abs(IntrpCumExctPtrnAll2-CumExtPtrnAll2(end)/2));
        StatJudge.ChCoG4AllDur2(nIntvl) = NchRes(Ncog2);       
        
    end;
    disp(sprintf('StatJudge.ChCoG4AllDur = %6.2f   %6.2f',StatJudge.ChCoG4AllDur(:)));
      
     %% Stat2: Vowel weighting
    CountVowelNQ= zeros(2,7);
    ListVowelNQ = {'a','i','u','e','o','N','Q'};
    for nIntvl = 1:NumIntvl % 2 interval 
      for nRpt = 1:SIMparam.NumRptWord,
        nwd = (nIntvl-1)*SIMparam.NumRptWord + nRpt;
        for nList = 1:length(ListVowelNQ)
            chr1 = char(ListVowelNQ(nList));
            CountVowelNQ(nIntvl,nList) = CountVowelNQ(nIntvl,nList) ...
                    + length(findstr(chr1,char(TextWord_SndPlay(nwd))));
        end;
      end;
    end;
    %CountVowelNQ
    if mean(sum(CountVowelNQ,2)) ~= 2*4,
        warning('Something wrong. Is it 4 morae word? or my mistake.');
    end;

    StatV = TableCoG_VowelMonoSyl();
    for nIntvl = 1:2
        CountVowel = sum(CountVowelNQ(nIntvl,1:5));
        StatJudge.ChCoG_StatVowelWeight(nIntvl) = ...
            (CountVowelNQ(nIntvl,1:5)*StatV.ChCoG_Vowel)/CountVowel;
    end;
    disp(sprintf('StatJudge.ChCoG_StatVowelWeight  = %6.2f   %6.2f',StatJudge.ChCoG_StatVowelWeight(:)));
    StatJudge.StatV = StatV;
    
   %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %   Stat6: Vowel CoG calculation %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  add 1 Nov 2014 (KY)
   
   StatJudge.ExctMeanCrossCorPtrn = zeros(5,100,2);
   Caltsec = 0.025; % ’†‰›‚©‚ç50ms
   
   if SIMparam.SwStatJudge == 6 | SIMparam.SwStatJudge == 7,
       % add 1 Nov 2014 (KY)
       
       ListVowel = {'a','i','u','e','o'};
       ListLgVowel = {'a:','i:','u:','e:','o:'};
       
       %find(strcmp(vowellabel_Intvl(:,:,2),ListVowel(1)) == 1)
       %pause;
      
       for nIntvl = 1:NumIntvl,  
         for nVowel = 1:5,
             NumMoraListVwl = find(strcmp(vowellabel(:,:,nIntvl),ListVowel(nVowel)) == 1 ...
                 | strcmp(vowellabel(:,:,nIntvl),ListLgVowel(nVowel)) == 1);
             NumMoraAllListVwl = find(strcmp(vowellabel,ListVowel(nVowel)) == 1 ...
                 | strcmp(vowellabel,ListLgVowel(nVowel)) == 1);
           
           if length(NumMoraAllListVwl) <= 1 | length(NumMoraListVwl) < 1,
             
             %disp('No-Vowel (or one-vowel) in this section -> No Calculation');
             StatJudge.ChCoG_Vwl(:,nVowel) = NaN;
             continue;
           else
             MtchVwl = NumMoraListVwl;
             if nIntvl == 2,
                 MtchVwl = MtchVwl + length(vowellabel(:,1,1))*2;
             end
             
             %% •ê‰¹‹æŠÔ‚Ì’²®
             if vowelStrtTsec(MtchVwl) > (vowelCTsec(MtchVwl) - Caltsec),
                 bigf = round(vowelStrtTsec(MtchVwl)/SGparam.Tshift);
             else
                 bigf = round((vowelCTsec(MtchVwl) - Caltsec)/SGparam.Tshift);
             end;
             
             if vowelEdTsec(MtchVwl) < (vowelCTsec(MtchVwl) + Caltsec),
                 endf = round(vowelEdTsec(MtchVwl)/SGparam.Tshift);
             else
                 endf = round((vowelCTsec(MtchVwl) + Caltsec)/SGparam.Tshift);
             end;
           
           end; % length(NumMoraAllListVwl) <= 1,
           
           ExctPtrn = zeros(NumCh,length(bigf));
           
           for nRpt6 = 1:length(bigf)
               ExctVwlsection = bigf(nRpt6):endf(nRpt6);
               ExctPtrn(:,nRpt6)= rms(GCSpecGram(:,ExctVwlsection,nIntvl),2);
           end;
           
           if length(bigf) > 1,
               ExctPtrnvwl = mean(ExctPtrn');
           elseif length(bigf) == 1
               ExctPtrnvwl = ExctPtrn;
           else
               error('No-Vowel (or one-vowel) in this section');
               ExctPtrnvwl = zeros(1,100);
           end;
           
           
           StatJudge.ExctMeanCrossCorPtrn(nVowel,:,nIntvl) = ExctPtrnvwl;
           
           CumExtPtrnvwl = cumsum(ExctPtrnvwl);
           ChIntrp = 0.01;
           NchRes = 1:ChIntrp:NumCh;
           IntrpCumExctPtrnvwl = interp1(1:NumCh,CumExtPtrnvwl,NchRes);
           [dummy Ncog3] = min(abs(IntrpCumExctPtrnvwl-CumExtPtrnvwl(end)/2));
           %StatJudge.ChCoG4VwlDur(nIntvl) = NchRes(Ncog3);
           StatJudge.ChCoG_Vwl(nIntvl,nVowel) = NchRes(Ncog3);
         end;
         
         disp(sprintf('StatJudge.ChCoG_Vwl = %6.2f  ',StatJudge.ChCoG_Vwl(nIntvl,:)));
         disp(' ')
         
         %disp(' ');
         %pause;
       end;
   end; %SeStatJudge = 6;
   
    %% Judgement based on Stat
    RspBtn = 0; % No response
    if     SIMparam.SwStatJudge == 0,
        ValJudge(1) = rand(1); % simple random for debug
        ValJudge(2) = rand(1);
    elseif SIMparam.SwStatJudge == 1,
        ValJudge(1) = StatJudge.ChCoG4AllDur(1);
        ValJudge(2) = StatJudge.ChCoG4AllDur(2);
    elseif SIMparam.SwStatJudge == 2,
        ValJudge(1) = (StatJudge.ChCoG4AllDur(1)-1)/(StatJudge.ChCoG_StatVowelWeight(1)-1);
        ValJudge(2) = (StatJudge.ChCoG4AllDur(2)-1)/(StatJudge.ChCoG_StatVowelWeight(2)-1);
    elseif SIMparam.SwStatJudge == 3,
        ValJudge(1) = StatJudge.ChCoG4AllDur(1) - StatJudge.ChCoG_StatVowelWeight(1);
        ValJudge(2) = StatJudge.ChCoG4AllDur(2) - StatJudge.ChCoG_StatVowelWeight(2);
    elseif SIMparam.SwStatJudge == 4, % several trial
        %ValJudge(1) = StatJudge.ChCoG4AllDur(1) + StatJudge.ChCoG_StatVowelWeight(1);
        %ValJudge(2) = StatJudge.ChCoG4AllDur(2) + StatJudge.ChCoG_StatVowelWeight(2); 
        %ValJudge(1) = StatJudge.ChCoG4AllDur2(1);  % NG
        %ValJudge(2) = StatJudge.ChCoG4AllDur2(2);  % NG
        %ValJudge(1) = StatJudge.ChCoG4AllDur2(1)- StatJudge.ChCoG_StatVowelWeight(1);  % NG
        %ValJudge(2) = StatJudge.ChCoG4AllDur2(2)- StatJudge.ChCoG_StatVowelWeight(2);  % NG        
        %ValJudge(1) = StatJudge.ChCoG4AllDur(1) - StatJudge.ChCoG_StatVowelWeight(1);
        %ValJudge(2) = StatJudge.ChCoG4AllDur(2) - StatJudge.ChCoG_StatVowelWeight(2)+20; 
        %ValJudge(1) = StatJudge.ChCoG4AllDur(1) - sqrt(StatJudge.ChCoG_StatVowelWeight(1));
        %ValJudge(2) = StatJudge.ChCoG4AllDur(2) - sqrt(StatJudge.ChCoG_StatVowelWeight(2)); 
        %VecLinModel = [ 156.36, -0.5535, 0.7839, 0.0197, -0.4187]/100; 
        VecLinModel = [12.7212   -1.1069    1.5678    0.0393   -0.8374]/100; 
        ValJudge(1) = 0;
        ValJudge(2) = VecLinModel*[1; StatJudge.ChCoG4AllDur(1);  StatJudge.ChCoG4AllDur(2); ...
                      StatJudge.ChCoG_StatVowelWeight(1); StatJudge.ChCoG_StatVowelWeight(2)];
        VecLinModelStep = [559.3251   -1.0573    1.6462  -13.9470  -14.8587    0.3546]/100;         
        ValJudge(2) = VecLinModelStep*[1; StatJudge.ChCoG4AllDur(1);  StatJudge.ChCoG4AllDur(2); ...
                      StatJudge.ChCoG_StatVowelWeight(1); StatJudge.ChCoG_StatVowelWeight(2); ...
                      StatJudge.ChCoG_StatVowelWeight(1)*StatJudge.ChCoG_StatVowelWeight(2); ];      
    elseif SIMparam.SwStatJudge == 5, % several trial
        ValJudge(1) = (StatJudge.ChCoG4AllDur(1) + StatJudge.ChCoG_StatVowelWeight(1))/StatJudge.ChCoG_StatVowelWeight(2);
        ValJudge(2) = (StatJudge.ChCoG4AllDur(2) + StatJudge.ChCoG_StatVowelWeight(2))/StatJudge.ChCoG_StatVowelWeight(1);
    elseif SIMparam.SwStatJudge == 6, % several trial
        
        ValNNan1 = find(StatJudge.ChCoG_Vwl(1,:)>0);
        ValNNan2 = find(StatJudge.ChCoG_Vwl(2,:)>0);
        
        %StatJudge.ChCoG4VwlDur(1) = prod(StatJudge.ChCoG_Vwl(1,ValNNan1));
        %StatJudge.ChCoG4VwlDur(2) = prod(StatJudge.ChCoG_Vwl(2,ValNNan2));
        
        StatJudge.ChCoG4VwlDur(1) = prod(StatJudge.ChCoG_Vwl(1,ValNNan1))/length(ValNNan1)
        StatJudge.ChCoG4VwlDur(2) = prod(StatJudge.ChCoG_Vwl(2,ValNNan2))/length(ValNNan2)
        
        %StatJudge
        ValJudge(1) = StatJudge.ChCoG4VwlDur(1);
        ValJudge(2) = StatJudge.ChCoG4VwlDur(2);
    elseif SIMparam.SwStatJudge == 7, %cross-corration
        
        StatJudge.ExctAllCrossCorPtrn
        Intvl1ClsCor = StatJudge.ExctMeanCrossCorPtrn(4,:,1);
        Intvl2ClsCor = StatJudge.ExctMeanCrossCorPtrn(4,:,2);
        
        [CC,lags] = xcorr(Intvl1ClsCor,Intvl2ClsCor,100,'coeff');
        %stem(lags(101:end),xc(101:end),'markerfacecolor',[0 0 1])
        lags(find(CC == max(abs(CC))))
        
        if lags(find(CC == max(CC)))>0,
            ValJudge(1) = 1;
            ValJudge(2) = 0;
        else
            ValJudge(1) = 0;
            ValJudge(2) = 1;
        end;
    
    else
        error('Not Prepared yet');
    end;
    RspBtn = 2 - (ValJudge(1) > ValJudge(2));   % Logic [1, 0] --> RspBtn [1, 2];
    
    disp(sprintf('ValJudge               = %6.2f   %6.2f',ValJudge(:)));
    StatJudge.ValJudge = ValJudge;
    
    disp(['Respone Button = ' int2str(RspBtn) ]);
    disp(['--------------------------']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Finalize     %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    RsltSim.StatJudge = StatJudge;
    RsltSim.RspBtn    = RspBtn;
    RsltSim.GCparam   = GCparam;
    RsltSim.SGparam   = SGparam;
    RsltSim.SIMparam  = SIMparam;
    
    %pause;
    
return;


%
%
function StatV = TableCoG_VowelMonoSyl();
    % Calculation by KY
    %
   StatV.ChCoG_MonoSyl = [ ...
   48.1400   48.4300   47.6100   47.5400   51.7700   51.7400   49.2500   55.1600   52.5600   45.9700;
   27.9000   28.6500   39.1900   30.5000   22.7200   26.8200   25.8700       NaN   27.8600       NaN;
   20.2700   21.8900   27.8500   28.3400   23.9400   22.2500   25.2900   24.6300   22.8700       NaN;
   61.6400   61.7800   57.4300   52.9500   59.8000   65.0700   64.5000       NaN   59.0400       NaN;
   36.2500   37.3800   43.9100   40.0300   38.4200   37.2600   38.4200   44.7300   43.7900       NaN;];
   StatV.ChCoG_Vowel = nanmean(StatV.ChCoG_MonoSyl,2);
   
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trash %%%%%%%%%%%%%%%%%%%%
         % [dummy AmpdB] = EqlzDigitalLevel(SndSim1,fs,ExdSndDigitalLeveldB,'LAeq');
        %RovingLeveldB = AmpdB(2);
        %SigSPLdB = ExpSndLeveldB + RovingLeveldB;
        %disp(['Signal Level = ' num2str(SigSPLdB) ', Roving Level = ' num2str(RovingLeveldB)]);
        %SIMparam.SigSPLdB(nIntvl) = SigSPLdB;

