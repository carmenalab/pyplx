function create_matfiles_from_plx_3(date,block)
% using Amy's plx2Mat function

%Date and block are strings
%dates={'071508','071608','071908'};
%dates = {'071508'};
%blocks  = {'a','a','a'};

subject = 'seba';
% e.g.
% plxF = 'C:\bmi_data\seba_lfp_data\seba042213\map_data\seba042213a.plx'
% saveDir = 'C:\bmi_data\seba_lfp_data\seba042213\map_data\';

datIncl = {'sig','wf','lfp','strobed'};

%plxF    = ['\\VBOXSVR\Matlab\map_data\' subject date block '.plx'];
plxF    = ['/work/pkhanna/convert_plx/plx_files/' subject date block '.plx'];
    
%saveDir = ['\\VBOXSVR\Matlab\map_data\'];
saveDir = ['/work/pkhanna/convert_plx/mat_files'];
        
matTag = '_all.mat';
        
plx2Mat(plxF, saveDir, matTag, datIncl);
        
disp(strcat('done with day ',date,', block ',block))
end

