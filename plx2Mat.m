function plx2Mat(plxF, saveDir, matTag, datIncl,inclSorted)

%plx2Mat(plxF, saveDir, matTag, datIncl)
%
%loads data from a plx file and saves it into a matlab file. 
%Will save a file in the specified path, with the name {plx file name}{matTag}
%
%input: plxF    - string with plx file name (full path)
%       saveDir - string with path for where to save the matlab file
%       matTag  - string with tag to append to matlab file name (default: '_raw.mat')
%       datIncl - string or cell of strings with names of data to include.
%                 possible datIncl:  sig     - spike times
%                                    wf      - spike wave-forms
%                                    kin     - kinematic channels
%                                    lfp     - lfp AD channels
%                                    emg     - emg AD channels
%                                    strobed - events
%	inclSorted - optional input string 'yes' will include unsorted units. 
%				If left empty or 'no', they will not be included,  
%note: if specified saveDir does not exist, function will create it. 
%
%requires the plexon OffineSDK functions available for download on plexon's website
%
% updated on 1-15-13, to include option to ignor unsorted units

kinChan  = (33:40);        
emgChan  = [50 52 56 61 63 64]; 
lfpChan  = (1:128)+64;
sigChan  = 1:128;


if exist('inclSorted','var')
	if strcmp(lower(inclSorted),'yes')
		units    = [1:4 0];  %sorted = 1:4, unsorted = 0
	elseif strcmp(lower(inclSorted),'no')
		units = [1:4];
	else
		error('Unrecognized inclSorted input')
	end
else
	units = [1:4];
end

unitTags = {'a', 'b', 'c', 'd', 'i'};
sigTag   = 'sig';
adTag    = 'AD';

if ~exist('matTag', 'var')
    matTag = '_raw.mat';
end
if ~exist('datIncl', 'var')
    datIncl = {'sig', 'kin', 'strobed'};
end
if ~iscell(datIncl)
    datIncl = {datIncl};
end
nDat = length(datIncl);


%check formatting of files/dirs
if ~strcmp(plxF(end-3:end), '.plx')
    plxF = strcat(plxF, '.plx');
end
if ~strcmp(saveDir(end), '/') && ~strcmp(saveDir(end), '\')
    saveDir = strcat(saveDir, '/');
end


%get file info about active spikes from plx header
addpath(genpath('/work/pkhanna/'));
% 
% disp(pwd)
% path

tscounts   = plx_info(plxF, 0);
activeChan = tscounts'>0; %shape to chan x unit
activeChan = activeChan(2:end,:); %tscounts shifted by 1
activeChan = circshift(activeChan, [0 -1]); %rotate unsorted to the last unit
saveVars   = cell(nDat,1);
cnt        = 1;
for i=1:nDat
    
    switch lower(datIncl{i})
        
        %get sig time-stamps
        case{'sig'}
            
            %loop through all possible channels/units
            for j=sigChan
                ch = num2str(j);
                ch = strcat(repmat('0', 1, 3-length(ch)), ch);
                for k=1:length(units)
                disp(strcat('getting ',num2str(k),'th unit'))                   
                    %only get data if active channel or don't know active channels yet
                    if activeChan(j,k)
                        [~, ts] = plx_ts(plxF, j,units(k));
                        
                        %only save if active spike data
                        if ~isequal(ts, -1)
                            eval([sigTag, ch, unitTags{k}, ' = ts;'])
                        end
                    end
                end
            end
            
            clear ts
            
            saveVars{cnt} = strcat(sigTag, '*');
            cnt           = cnt+1;
            
        case{'wf'}
            
            %loop through all possible channels/units
            for j=sigChan
                ch = num2str(j);
                ch = strcat(repmat('0', 1, 3-length(ch)), ch);
                for k=1:length(units)
                    disp(strcat('getting', num2str(k),'th waveform'))
                    %only get data if active channel or don't know active channels yet
                    if activeChan(j,k)
                        [~, ~, ts, wave] = plx_waves_v(plxF, j,units(k)); %#ok<NASGU>
                        
                        %only save if active spike data
                        if ~isequal(ts, -1)
                            eval([sigTag, ch, unitTags{k}, '_wf = wave;'])
                            eval([sigTag, ch, unitTags{k}, '_wf_ts = ts;'])
                        end
                    end
                end
            end
            clear wave ts
            
            saveVars{cnt} = strcat(sigTag, '*');
            cnt           = cnt+1;
        
        case{'kin'}
            
            %loop through kin ad channels
            for j=kinChan
                [~, ~, ts, ~, ad] = plx_ad_v(plxF, j-1);   %#ok<NASGU> %ad channel indexing shifted by 1
                
                if ~isequal(ts, -1)
                    eval([adTag, num2str(j), ' = ad;'])
                    eval([adTag, num2str(j), '_ts = ts;'])
                end
            end
            clear ad ts
            
            saveVars{cnt} = strcat(adTag, '*');
            cnt           = cnt+1;
            
        case{'lfp'}
            
            %loop through lfp chanels
            for j=lfpChan
                [~, ~, ts, ~, ad] = plx_ad_v(plxF, j-1);   %#ok<NASGU> %ad channel indexing shifted by 1
                
                if ~isequal(ts, -1)
                    eval([adTag, num2str(j), ' = ad;'])
                    eval([adTag, num2str(j), '_ts = ts;'])
                end
            end
            clear ad ts
            
            saveVars{cnt} = strcat(adTag, '*');
            cnt           = cnt+1;
            disp('done with lfp');
            
        case{'emg'}
            
            %loop through emg chanels
            for j=emgChan
                [~, ~, ts, ~, ad] = plx_ad_v(plxF, j-1);   %#ok<NASGU> %ad channel indexing shifted by 1
                
                if ~isequal(ts, -1)
                    eval([adTag, num2str(j), ' = ad;'])
                    eval([adTag, num2str(j), '_ts = ts;'])
                end
            end
            clear ad ts
            
            saveVars{cnt} = strcat(adTag, '*');
            cnt           = cnt+1;
            
        case{'strobed'}
            
            [~, ts, sv] = plx_event_ts(plxF, 257);
            if ~isequal(ts, -1)
                Strobed = [ts sv]; %#ok<NASGU>
            end
            
            saveVars{cnt} = 'Strobed';
            cnt           = cnt+1;
            disp('done with strobed');
            
        otherwise
            warning(['data type not recognized: ', datIncl{i}]) %#ok<WNTAG>
            
    end
    
end %data types


%get file-name for mat file
dash = cat(2, strfind(plxF, '\'), strfind(plxF, '/'));
dash = sort(dash);
if isempty(dash)
    start = 1;
else
    start = dash(end)+1;
end
fbase = plxF(start:end-4);
disp('writing, part1')

if ~exist(saveDir, 'dir')
    mkdir(saveDir)
    fprintf('Creating directory: %s\n', saveDir)
end
fprintf('Saving %s\n', strcat(saveDir, fbase, matTag))

save(strcat(saveDir, fbase, matTag), saveVars{:})

disp('done saving')

    
            
                    
            
            
