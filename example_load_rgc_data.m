
function expdata = example_load_rgc_data(varargin)
%
%
%%% example_load_rgc_data %%%
%
%
% This function loads the data for any example cell from the accompanying
% dataset of Khani and Gollisch (2021). The dataset is available at:
% https://gin.g-node.org/gollischlab/Khani_and_Gollisch_2021_RGC_spike_trains_chromatic_integration
% This function first get the path of the downloaded database and then loads
% the frametimes and spike times of a selected ganglion cells from the list
% of avaiable good cells. The list that can be found at list_of_good_cells.txt
% inside each folder of the database.
% example use: 1- example_load_rgc_data();
%              2- example_load_rgc_data('datapath','PATH TO DOWNLOAD FOLDER')
%              3- example_load_rgc_data('datapath','PATH TO DOWNLOAD FOLDER', stimulusnumber, 7)
%              4- example_load_rgc_data('datapath','PATH TO DOWNLOAD FOLDER', cellnumber, 12, clusternumber, 12)
%              4- example_load_rgc_data('datapath','PATH TO DOWNLOAD FOLDER', stimulusnumber, 7, cellnumber, 12, clusternumber, 12)
%
%
% ===============================Inputs====================================
%
%   datapath        :   path to the folder in the database.
%   stimulusnumber  :   number to select the correct stimulus.
%   cellnumber      :   an example cell number from the list of good cells.
%   clusternumber   :   a cluster number for the selected cell.
%
%================================Output====================================
%
%   expdata         :   spiketime, frametimes, stimulus parameters and cell
%                       id information of the selected cell.
%
% written by Mohammad, 18.02.2021

%--------------------------------------------------------------------------------------------------%
%----------                              main-function                                   ----------%
%--------------------------------------------------------------------------------------------------%

% first load the data for an example cell
[dp, stimulusnumber, cellnumber, clusternumber]         =       select_example_cell(varargin{:});
% now load the date that selected example cell
expdata         =       load_data_for_selected_cell(dp, stimulusnumber, cellnumber, clusternumber);
expdata.datapath = dp;
expdata.stimulusnumber = stimulusnumber;
expdata.cellnumber = cellnumber;
expdata.clusternumber = clusternumber;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [dp, stimulusnumber, cellnumber, clusternumber] = select_example_cell(varargin)
% first parse the user inputs
p = inputParser();      % check the user options.
p.addParameter('datapath', [],  @(x) ischar(x) || isstring(x));
p.addParameter('stimulusnumber', [], @isnumeric);
p.addParameter('cellnumber', [], @isnumeric);
p.addParameter('clusternumber', [], @isnumeric);
p.addParameter('help', false,  @(x) islogical(x) || (isnumeric(x) && ismember(x,[0,1])));

p.parse(varargin{:});
% defualt parameters
defpara     =    p.Results;

if isempty(defpara.datapath)
    dp                  =       uigetdir('Select an experiment folder:');
else
    dp                  =       defpara.datapath;
end

if isempty(defpara.cellnumber) || isempty(defpara.clusternumber)
    listofgoodcells     =       load([dp,filesep,'list_of_good_cells.txt'],'-ascii');
    cellabels   =   sprintfc('cell %g, cluster %g',listofgoodcells);
    if numel(cellabels) > 1
        idx         =   listdlg('PromptString',{'Select an example ganglion cell:'},...
            'SelectionMode','single','ListString',cellabels,'listsize',[250 250]);
    else
        idx         =   1;
    end
    cellnumber          =       listofgoodcells(idx,1);
    clusternumber       =       listofgoodcells(idx,2);
else
    cellnumber      =   defpara.cellnumber;
    clusternumber   =   defpara.clusternumber;
end

if isempty(defpara.stimulusnumber)
    stimnames           =       importdata([dp,filesep,'stimuli_names.txt']);
    cinames     =   stimnames(~(cellfun('isempty',(strfind(stimnames,'_Chromatic_Integration_')))));
    cinames     =   cinames((cellfun('isempty',(strfind(cinames,'Local')))));
    cinames     =   cinames((cellfun('isempty',(strfind(cinames,'Grating')))));
    
    if numel(cinames) > 1
        idx     = listdlg('PromptString',{'Select an example dataset:'},...
            'SelectionMode','single','ListString',cinames,'listsize',[450 100]);
        cinames     =   cinames{idx};
    end
    stimulusnumber      =       str2double(extractBefore(cinames,'_Chromatic'));
else
    stimulusnumber  =       defpara.stimulusnumber;
end

end

%--------------------------------------------------------------------------------------------------%

function expdata = load_data_for_selected_cell(dp, stimnum, cellnum, clusternum)
%
ftfolder        =       [dp,filesep,'frametimes'];
if ~exist(ftfolder,'dir'),      ftfolder    =   dp;     end
rasfolder       =       [dp,filesep,'spiketimes'];
if ~exist(rasfolder,'dir'),      rasfolder    =   dp;     end
paramfolder     =       [dp,filesep,'stimulusparameters'];
if ~exist(paramfolder,'dir'),      paramfolder    =   dp;     end

ftnames         =       dir([ftfolder, filesep, num2str(stimnum,'%02g'),'*_frametimings.txt']);
rastername      =       [rasfolder, filesep, sprintf('%d_SP_C%d%02d.txt', stimnum, cellnum, clusternum)];

warning('off','all');   % this is supress the warning about long names, hopefully nothing happen in between these lines!
expdata.spiketimes      =       load(rastername, '-ascii');
expdata.frametimes      =       load([ftfolder,filesep,ftnames.name], '-ascii');

stimpath                =       dir([paramfolder,filesep,num2str(stimnum,'%02g'),'*_parameters.txt']);

% now reading the txt file with parameters
fid             =       fopen([paramfolder,filesep,stimpath.name]);
tline           =       fgetl(fid);
while ischar(tline)
    tline       =       fgetl(fid);
    if tline    ==      -1,     break;      end
    fn          =       extractBefore(tline,' = ');
    if isempty(fn),     continue;           end
    val         =       extractAfter(tline,' = ');
    % to convert to double
    if ~isnan(str2double(val))
        val     =       str2double(val);
    end
    % to convert to logical
    if strcmp(val,'true'),        val = true;           end
    if strcmp(val,'false'),       val = false;          end
    % store the values in a structure
    stimpara.(fn)       =       val;
end
fclose(fid);

expdata.stimpara    =       stimpara;
expdata.name        =       extractBefore(ftnames.name,'_frametimings.txt');
expdata.paraname    =       extractBefore(stimpath.name,'_parameters.txt');
if strcmp(expdata.name , expdata.paraname)
    expdata         =       rmfield(expdata,'paraname');
end
warning('on','all');

end

%--------------------------------------------------------------------------------------------------%

