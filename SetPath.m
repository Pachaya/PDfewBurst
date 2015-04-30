function [ PATH ] = SetPath(  )
%Get the computername and set the path
%   VSLAB-PC = Pachaya's machine at CMS#307
%   USER-PC = VS lab's machine at experimental room (CMS#306)

[s,computername] = dos('ECHO %COMPUTERNAME%');
computername = computername(1:end-1); % Remove '\n' at the end
% VSLAB = sprintf('VSLAB-PC\n');
% EXP_PC = sprintf('USER-PC\n');
disp('--------------------------------------------------------------------')
disp(['Computer name :: ' computername ])
switch computername
        case 'VSLAB-PC'
           PATH = 'E:\PDmodelFewBurst\SimResult\';
        case 'USER-PC'
           PATH = 'D:\Pachaya\150313 Code\SimResult\'; 
    otherwise
            disp('Unknown Computer name :: Use current directory for path')
            disp(['Current Directory : ' pwd]);
            PATH = '';
            
end
   disp(['Setting Path to : ' PATH ]);
disp('--------------------------------------------------------------------')

end

