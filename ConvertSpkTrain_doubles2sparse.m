%% Convert spktrain from doubles to sparse
clear all

loc = ['E:\PDmodelFewBurst\SimResult\OscInput_varyTCtype_Sim\CompareConTypes_50Hz__rTC_50_250_wmTC_50_500_OSC_F_20_40_OSC_amp_0_1_Wspk_0.001_0.001\'];
fname = ['Activity_avgPnegecpW_27-May-2015.mat'];
load([loc fname]);
fname2 = 'Activity_avgPnegexpW';
s = size(ACT_Record);

for p1 = 1: s(1)
    for p2  = 1: s(2)
        for p3  = 1: s(3)
            for p4  = 1: s(4)
                for p5  = 1: s(5)
                    fprintf('------------------------------\n')
                    fprintf('------ %d %d %d %d %d --------\n',p1,p2,p3,p4,p5 )
                    fprintf('------------------------------\n')
                    
                    layers = fieldnames(ACT_Record{p1,p2,p3,p4,p5});
                    for ll = 1 : length(layers)
                        fprintf('\t%s\n', layers{ll})
                        CaT = fieldnames(ACT_Record{p1,p2,p3,p4,p5}.(layers{ll}));
                        for cc = 1 : length(CaT)
                            fprintf('\t\t%s\n',CaT{cc});
                            ctypes = fieldnames(ACT_Record{p1,p2,p3,p4,p5}.(layers{ll}).(CaT{cc}));
                            for tt = 1 : length(ctypes) -1  % E, I ,All
                                ACT_Record{p1,p2,p3,p4,p5}.(layers{ll}).(CaT{cc}).(ctypes{tt}).spktrain = sparse( ACT_Record{p1,p2,p3,p4,p5}.(layers{1}).(CaT{cc}).(ctypes{tt}).spktrain );
                                fprintf('\t\t\t%s\n', ctypes{tt})
                            end
                        end
                    end
                    
                    
                end
            end
        end
    end
end

%%

save([loc 'sparse' fname2  date '.mat' ], 'ACT_Record','PARAMETERS','-v7.3');
