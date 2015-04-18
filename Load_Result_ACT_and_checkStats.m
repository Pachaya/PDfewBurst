% Download Activity Result and make statistical Test'


%%
% load('Result_ACT_Record_PDfewBurst_GPmVLmd1.mat')


ACT_Record_Test = cell(TRIAL_NO,1);
ncells = 1150;
for TRIAL_NO =  1: NUM_TRIAL
    % T - test
    Basal_Act = ACT_Record{TRIAL_NO};
    WT = Basal_Act.WT; KO = Basal_Act.KO;
    disp('--------------------------------------------------------------------------------------------------')
    disp([' Trial # = ' num2str(TRIAL_NO)]);
    disp('--------------------------------------------------------------------------------------------------')
    disp('E cells : two-tailed t-test')
    [hE,pE] = ttest2(WT.All.fr_data, KO.All.fr_data, 'Tail','both');
    disp([ ' p-value = ' num2str(pE)])
    if(hE)
        disp('H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)')
    else
        disp('Do not rejected H0 : average normal firing rate of E cells in WT does not significantly different from KO (p > 0.05)')
    end
    Basal_Act.ttest2.hE = hE;
    Basal_Act.ttest2.pE = pE;
    
    
    % Sample T-Test 
      
        Nsample = 15;
        Nrepeat = 100;
        disp(['===== Sample Test : #sample = ' num2str(Nsample) ' for ' num2str(Nrepeat) ' trials ====='])
        tresult = cell(Nrepeat,1);
        plist = zeros(Nrepeat,1);
        hlist = zeros(Nrepeat,1);
        tresult_t = cell(Nrepeat,1);
        plist_t = zeros(Nrepeat,1);
        hlist_t = zeros(Nrepeat,1);
        
        for ii = 1: Nrepeat
            sample = randperm(ncells, Nsample);
            [ht,pt] = ttest2(WT.All.fr_data(sample), KO.All.fr_data(sample));
            [h,p] = ranksum(WT.All.fr_data(sample), KO.All.fr_data(sample));
            
            tresult{ii}.sampleID = sample;
            tresult{ii}.h = h;
            tresult{ii}.p = p;
            %     disp(['h:' num2str(h) ', p:' num2str(p)])
            plist(ii) = p;
            hlist(ii) = h;
            
            tresult_t{ii}.sampleID = sample;
            tresult_t{ii}.h = ht;
            tresult_t{ii}.p = pt;
            %     disp(['ht:' num2str(ht) ', p:' num2str(pt)])
            plist_t(ii) = pt;
            hlist_t(ii) = ht;
        end
        
        disp('##Ranksum##')
        %  sum((plist>0.05))
        disp([num2str(sum((plist>0.05))) ' cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)'])
        disp('##T-test#')
        %  sum((plist_t>0.05))
        disp([num2str(sum((plist_t>0.05))) ' cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)'])
        % Common cases
        disp('The other cases are rejecting H0 : average normal firing rate of all cells in WT significantly different from KO (p < 0.05)')
        disp('')
        disp('common cases in Ranksum and T-test')
        find((plist>0.05) == (plist_t>0.05))
        Basal_Act.sampleTest.ranksum = tresult;
        Basal_Act.sampleTest.ttest = tresult_t;
        ACT_Record_Test{TRIAL_NO} = Basal_Act;
end




%% 
%{
--------------------------------------------------------------------------------------------------
 Trial # = 1
--------------------------------------------------------------------------------------------------
E cells : two-tailed t-test
 p-value = 2.3088e-14
H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)
===== Sample Test : #sample = 15 for 100 trials =====
##Ranksum##
10 cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)
##T-test#
87 cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)
The other cases are rejecting H0 : average normal firing rate of all cells in WT significantly different from KO (p < 0.05)
common cases in Ranksum and T-test

ans =

     6
     9
    27
    35
    39
    40
    45
    52
    60
    62
   100

--------------------------------------------------------------------------------------------------
 Trial # = 2
--------------------------------------------------------------------------------------------------
E cells : two-tailed t-test
 p-value = 1.282e-12
H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)
===== Sample Test : #sample = 15 for 100 trials =====
##Ranksum##
11 cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)
##T-test#
90 cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)
The other cases are rejecting H0 : average normal firing rate of all cells in WT significantly different from KO (p < 0.05)
common cases in Ranksum and T-test

ans =

    20
    31
    37
    51
    61
    88
    94

--------------------------------------------------------------------------------------------------
 Trial # = 3
--------------------------------------------------------------------------------------------------
E cells : two-tailed t-test
 p-value = 4.0231e-15
H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)
===== Sample Test : #sample = 15 for 100 trials =====
##Ranksum##
10 cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)
##T-test#
87 cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)
The other cases are rejecting H0 : average normal firing rate of all cells in WT significantly different from KO (p < 0.05)
common cases in Ranksum and T-test

ans =

    17
    20
    36
    47
    56
    60
    65
    85
    90

--------------------------------------------------------------------------------------------------
 Trial # = 4
--------------------------------------------------------------------------------------------------
E cells : two-tailed t-test
 p-value = 1.8758e-11
H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)
===== Sample Test : #sample = 15 for 100 trials =====
##Ranksum##
13 cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)
##T-test#
91 cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)
The other cases are rejecting H0 : average normal firing rate of all cells in WT significantly different from KO (p < 0.05)
common cases in Ranksum and T-test

ans =

     8
    41
    60
    75

--------------------------------------------------------------------------------------------------
 Trial # = 5
--------------------------------------------------------------------------------------------------
E cells : two-tailed t-test
 p-value = 7.414e-13
H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)
===== Sample Test : #sample = 15 for 100 trials =====
##Ranksum##
7 cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)
##T-test#
90 cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)
The other cases are rejecting H0 : average normal firing rate of all cells in WT significantly different from KO (p < 0.05)
common cases in Ranksum and T-test

ans =

    30
    44
    58
    68
    84
    
    }%