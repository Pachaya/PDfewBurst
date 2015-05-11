function [ttest_res,ranksum_res, sum_sigDiffRS, sum_sigDiffTT] = sampleStatsTest(WTfr, KOfr, Nsample, Nrepeat)
% Sample Ranksum and T-Test for WT and KO 
%   WTfr, KOfr = average firing rate of WT cells and KO cells, respectively
%   Nsampl = number of sample
%   Nrepeat = total number of test
                              
%                         WTfr = Basal_Act.VL.WT.All.fr_data;
%                         KOfr = Basal_Act.VL.KO.All.fr_data;
%                         
%                         Nsample = 15;
%                         Nrepeat = 100;
                        ncells = min(length(WTfr),  length(KOfr));
                        disp(['===== Statistical Sample Test : #sample = ' num2str(Nsample) ' for ' num2str(Nrepeat) ' trials ====='])
                        rs_result = cell(Nrepeat,1);
                        plist = zeros(Nrepeat,1);
                        hlist = zeros(Nrepeat,1);
                        tt_result = cell(Nrepeat,1);
                        plist_t = zeros(Nrepeat,1);
                        hlist_t = zeros(Nrepeat,1);
                        
                        for ii = 1: Nrepeat
                            sample = randperm(ncells, Nsample);
                            [ht,pt] = ttest2(WTfr(sample), KOfr(sample));
                            [h,p] = ranksum(WTfr(sample), KOfr(sample));
                            
                            rs_result{ii}.sampleID = sample;
                            rs_result{ii}.h = h;
                            rs_result{ii}.p = p;
                            %     disp(['h:' num2str(h) ', p:' num2str(p)])
                            plist(ii) = p;
                            hlist(ii) = h;
                            
                            tt_result{ii}.sampleID = sample;
                            tt_result{ii}.h = ht;
                            tt_result{ii}.p = pt;
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
                        disp('-----------------------------------------------------------------------------------------------------------------------------')
                        
                        
             if nargout > 0
                 ranksum_res = rs_result;
                 ttest_res = tt_result;
                 sum_sigDiffRS = sum((plist>0.05));
                 sum_sigDiffTT = sum((plist_t>0.05));
             end      

end

