function [found] = CheckFileExist( dirLoc, Name_postfix  )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%   tmpFileExists{cell_type} = exist([dirLoc 'SomaVolt_' Name_postfix '.txt'], 'file');
                    if exist([dirLoc 'SomaVolt_' Name_postfix '.txt'], 'file')
                        disp(' ----- Yes -----')
                        %                         fname = [dirLoc 'SomaVolt_' Name_postfix '.txt'];
                        %                         aa = importdata(fname);
                        %                         size(aa)
                        found =1;
                    else
                        disp('### Sorry ###')
                        found = 0;
                    end
                    
end

