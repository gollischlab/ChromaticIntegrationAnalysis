

function Analyze_Chromatic_Integration_Stimulus(varargin)
%
%
%%% Analyze_Chromatic_Integration_Stimulus %%%
%
%
% This function analyzes the response of the ganglion cells to the chromatic
% integration stimulus. This is a supplement code for the Khani & Gollisch
% (2021) study of chromatic integration in the mouse retina.
% This function first loads the data in a format that is provided in by the
% supplementary database of the study. It then generate rasters and PSTHs for
% all the 22 contrast combinations that was shown to each retina. The order
% of the contrasts were randomized and it is recreated using Fisher-Yates
% random permutation algorithem. Check fisher_Yates_shuffle_order and
% reorder_contrast_spikes functions below for more information.
% The contrast values of both colors range from -20 to 20 in steps of 2% change.
% They are categorized in two groups of Green-On-UV-Off and Green-Off-UV-On.
%-------------------------------Green-On-UV-Off-----------------------------
% Green =>  20    18    16    14    12    10     8     6     4     2     0
% UV    =>   0    -2    -4    -6    -8   -10   -12   -14   -16   -18   -20
%-------------------------------Green-Off-UV-On-----------------------------
% Green => -20   -18   -16   -14   -12   -10    -8    -6    -4    -2     0
% UV    =>   0     2     4     6     8    10    12    14    16    18    20
%
%
% ===============================Inputs====================================
%
%   datapath        :   path to the folder in the database.
%   stimulusnumber  :   number to select correct stimulus.
%   cellnumber      :   an example cell number from the list of good cells.
%   clusternumber   :   a cluster number for the selected cell.
%
%================================Output====================================
%
%   chromatic integration curves and PSTHs for the selected cell.
%
% written by Mohammad, 12.02.2021

%--------------------------------------------------------------------------------------------------%
%----------                              main-function                                   ----------%
%--------------------------------------------------------------------------------------------------%

% first load the data for an example cell
[dp, stimulusnumber, cellnumber, clusternumber] = select_example_cell(varargin{:});
% get the parameters for the selected example cell
expdata     =       ci_parameters(dp, stimulusnumber, cellnumber, clusternumber);
% calculate the rasters and PSTH for the stimulus and the gray frame before the stimulus
ci_data     =       ci_analysis(expdata);
% plot chromatic integration curves and all the psths
plot_chromatic_integration(ci_data, cellnumber, clusternumber, expdata.stimpara);

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function expdata = ci_parameters(dp, stimulusnum, cellnum, clusternum,varargin)
%
expdata         =       load_data_for_selected_cell(dp, stimulusnum, cellnum, clusternum);
para            =       expdata.stimpara;
% loading clusters and frametimings
ft              =       expdata.frametimes;    % work in seconds

% getting contrasts of the experiment (-20:2:20)
contgr  =   round([para.mincontrast:para.contrastdiff:0,0:para.contrastdiff:para.maxcontrast]*100);
contuv  =   round([0:para.contrastdiff:para.maxcontrast,para.mincontrast:para.contrastdiff:0]*100);
% re-ordering stimulus for plotting and presentation purpose, here we start
% with on green and reverse off green instead of off to on green as in stimulus program.
para.contrastIdx        =       [length(contgr) : -1 : 1 + length(contgr)/2, 1 : length(contuv)/2];
para.contrastgreen      =       contgr(para.contrastIdx);
para.contrastuv         =       contuv(para.contrastIdx);
% generating contrast values
para.numtrials          =       floor( length( ft(2 : 2 : end)) / length(contgr) );
seed                    =       para.seed;
stimorder               =       nan(para.numtrials,length(contgr));

% getting the order and location of each contrast from the psudo random sequence
for ii      =   1 : para.numtrials
    [stimorder(ii,:),seed,~]    =   fisher_Yates_shuffle_order(seed,contgr);
end

% experiment parameters
refreshrate             =       60;
para.binlength          =       10 / 1e3;         % 10 ms each bin
para.delay              =       0;                % 0 ms delay
para.activeregion       =       [0.025, 0.25];    % region used to measure response (50 to 250 ms)
para.refreshrate        =       refreshrate;
% make a bin vector for stimulus, used for making PSTH
para.stimduration       =       (para.stimduration / refreshrate); % in sec
stimBin                 =       linspace(0, para.stimduration, (para.stimduration / para.binlength)+1);   % +1 here to avoid issues with histc last bin.
% make a bin vector for preframe, used for making PSTH
para.pfrduration        =       (para.preframes / refreshrate); % in sec
pfrBin                  =       linspace(0, para.pfrduration,(para.pfrduration / para.binlength) + 1);
% seting output
expdata.stimorder               =       stimorder;
expdata.stimulusbinvector       =       stimBin;
expdata.preframebinvector       =       pfrBin;
expdata.stimpara                =       para;

end

%--------------------------------------------------------------------------------------------------%

function ci_data = ci_analysis(expdata)
%
ft          =       expdata.frametimes;
spks        =       expdata.spiketimes;
para        =       expdata.stimpara;
% for whole experiment generally has 2000 ms preframe and 500 ms stimulus
[~, allpreframespikes] = spikes_per_frametimes(ft, spks, para.stimduration, para.pfrduration,para);
% this is to get rasters and psth for the graybackground before the stimulus
[ci_data.preframe.allpfrRasters, ci_data.preframe.allpfrPSTH] = reorder_contrast_spikes(allpreframespikes,...
    expdata.stimorder,expdata.preframebinvector,para);

%%% for 500 ms before and after stimulus
[stimulusspikes, preframespikes] = spikes_per_frametimes(ft,spks,para.stimduration,para.stimduration,para);
% rasters and pstrh for the 500 ms stimulus
[ci_data.stimulusRasters, ci_data.stimulusPSTH] = reorder_contrast_spikes(stimulusspikes,...
    expdata.stimorder,expdata.stimulusbinvector,para);
% rasters and pstrh for the 500 ms gray background before the onset of the stimulus
[ci_data.preframe.preframeRasters, ci_data.preframe.preframePSTH] = reorder_contrast_spikes(preframespikes,...
    expdata.stimorder,expdata.stimulusbinvector,para);
ci_data.preframePSTH = ci_data.preframe.preframePSTH;

end

%--------------------------------------------------------------------------------------------------%

function varargout = fisher_Yates_shuffle_order(seed,inputVec,varargin)
%
%%% fisher_Yates_shuffle_order %%%
%
%
% This function generate psudorandom permution similar to randperm in MATLAB
% but works with psudorandom number generator ran1. It also gives back the
% most recent seed value to continue the permuation in case of repeated trials.
% note that the direction of permutation is along x-axix or for columns of
% MATALB not for the rows.
%
%
% ===============================Inputs====================================
%
%   seed : seed value for random number generation.
%   inputVec : input vector used for permutation.
%
%================================Output====================================
%
%   testOrder : vector of permuted indices for the inputVec.
%   newSeed : the recent seed that used in the ran1 function.
%   outputVec : the permuted input vector along x-axis
%
% Note that the permution algorithem is based on Fisher-Yates shuffle
% algorithem identical to what is used in the stimulus program.
% for more info check :
% https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
%
% written by Mohammad, 01.02.2016

newSeed         =   seed;
testOrder       =   zeros(1,size(inputVec,2));
testOrder(1)    =   1;

for i   =   2:length(inputVec)-1
    
    [randVal,newSeed] = ran1(newSeed);
    
    j               =   ceil(i*randVal);    % based on Fischer-Yates algorithem
    testOrder(i)    =   testOrder(j);
    testOrder(j)    =   i;
    
end

testOrder           =   [testOrder(end),testOrder(1:end-1)]+1;   % to match MATLAB indexing
varargout{1}        =   testOrder;
varargout{2}        =   newSeed;
for j   =   1:size(inputVec,1)
    varargout{3}(j,:)   =   inputVec(j,testOrder);
end

end

%--------------------------------------------------------------------------------------------------%

function [stimSpk, pfrSpk] = spikes_per_frametimes(ft,spk,stimdur,pfrdur,para,varargin)
%
delay               =           para.delay / 1e3;    % change delay time to second to match ft and spk
[stimSpk,pfrSpk]    =           deal( cell( ceil( numel( ft )/ 2), 1));
iter                =           1;

for ii   =   2:2:numel(ft)
    
    stimSpk{iter}   = spk(and(spk > ft(ii)+delay,spk <= ft(ii)+(stimdur+delay))) - (ft(ii)+delay);
    % before stimulus onset to get preframe
    pfrSpk{iter}    = spk(and(spk > ft(ii)-(pfrdur+delay),spk <= ft(ii))) - (ft(ii)-(pfrdur+delay));
    
    if isempty(stimSpk{iter}),  stimSpk{iter}    =       NaN;        end
    if isempty(pfrSpk{iter}),   pfrSpk{iter}     =       NaN;        end
    iter            =           iter+1;
end

end

%--------------------------------------------------------------------------------------------------%

function varargout = reorder_contrast_spikes(spk,order,binVec,para,varargin)
%
neworder        =       reshape(order',1,[]);       % reshape stimulus order
spk             =       spk(1:length(neworder))';   % cutting the last unfinished stimulus series

% pre-allocating memory
spkorder        =       cell(size(order));
contspk         =       cell(1,max(neworder));
psth            =       nan(length(binVec),max(neworder));

for ii    =  1 : max(neworder)
    
    thisorder       =   neworder == ii;       % indexing the orders from 1 to end
    spkorder(:,ii)  =   spk( thisorder );     % sorting spikes based on previous order
    contspk{ii}     =   CelltoMatUE( spkorder(:,ii));
    if size(contspk{ii},2)  <  2
        psth(:,ii)  =   (histc(contspk{ii}',binVec)/size(contspk{ii},1))*(1/para.binlength);  %#ok to fix situation with 1 spike
    else
        psth(:,ii)  =   mean(histc(contspk{ii}',binVec),2)*(1/para.binlength);   %#ok must transpose contspk for histc to work properly.
    end
end
psth                =   psth(1:length(binVec)-1,:);  % delete last zero bin.
% Here we reorder the psth and spikes in a way that the opposing stimuli
% face each other. the indices 22:-1:11 are set to 1:11, these are green on
% to uv off. The indices 1:1:11 are set to 12 to 22, these are green off uv
% on. for the rest of the function no index correction is made.
% the idices for -20:2:20 are as follows
% for first 11 or Green-ON, UV-OFF
% Green =>  20    18    16    14    12    10     8     6     4     2     0
% UV    =>   0    -2    -4    -6    -8   -10   -12   -14   -16   -18   -20
% for second 11 or UV-ON, Green-OFF
% Green => -20   -18   -16   -14   -12   -10    -8    -6    -4    -2     0
% UV    =>   0     2     4     6     8    10    12    14    16    18    20

contidx             =   [length(para.contrastgreen) : -1:1 + length(para.contrastgreen)/2, 1: length(para.contrastuv)/2];

varargout{1}        =       contspk(contidx);
varargout{2}        =       psth(:,contidx);
varargout{3}        =       spkorder(:,contidx);

end

%--------------------------------------------------------------------------------------------------%

function [ci, ciang] = chromatic_integration_curves(stimpsth, pfrpsth, para, varargin)
%
contnum         =   size(stimpsth, 2) / 2;

timevec         =   linspace(0, para.stimduration, para.stimduration / para.binlength);
activtime       =   timevec >= para.activeregion(1) & timevec <= para.activeregion(2);

% only between 50-250 ms is counted!
spkdiff         =   stimpsth(activtime,:) - pfrpsth(1+end-sum(activtime):end,:);
respdiff        =   mean(spkdiff,1);    % difference in response
respsem         =   nan_sem(spkdiff);   % standard error of mean

grONuvOFF       =   respdiff(1 : contnum);
grOFFuvON       =   respdiff(1 + contnum : contnum*2);

grONuvOFFsem    =   respsem(1 : contnum);
grOFFuvONsem    =   respsem(1 + contnum : contnum*2);

gNuFspk         =   spkdiff(:, 1:contnum);
gFuNspk         =   spkdiff(:, 1+contnum : contnum*2);

ci.grONuvOFF            =       grONuvOFF;
ci.grONuvOFFsem         =       grONuvOFFsem;
ci.grONuvOFFspkdiff     =       gNuFspk;
ci.grOFFuvON            =       grOFFuvON;
ci.grOFFuvONsem         =       grOFFuvONsem;
ci.grOFFuvONspkdiff     =       gFuNspk;

if nargout > 1
    gNuFangles  =   [para.contrastgreen(1:contnum);para.contrastuv(1:contnum)];
    gFuNangles  =   [para.contrastgreen(1 + contnum:contnum*2); para.contrastuv(1 + contnum:contnum*2)];
    ciang.gNuFangles    =       gNuFangles;
    ciang.gFuNangles    =       gFuNangles;
end

end

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

function y = nan_sem(x,dim)
%
if isempty(x),      y = NaN;	return;     end
if nargin < 2,      dim = find(size(x)~=1, 1 );	if isempty(dim),    dim = 1; 	end; end
% Find NaNs in x and nanmean(x)
nans        =       isnan(x);
count       =       size(x,dim) - sum(nans,dim);
% Protect against a  all NaNs in one dimension
i           =       find(count==0);
count(i)    =       1;
y           =       nanstd(x,[],dim)./sqrt(count);
y(i)        =       i + NaN;

end

%--------------------------------------------------------------------------------------------------%

function OutputMat = CelltoMatUE (input_cell)
%
inptsize        =       size(input_cell);
maxLength       =       max(cellfun('length',input_cell));
if inptsize(2)  >       inptsize(1)
    OutputMat   = cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),input_cell,'un',0));
else
    OutputMat   = cell2mat(cellfun(@(x)cat(2,x',nan(1,maxLength-length(x))),input_cell,'un',0));
end

end

%--------------------------------------------------------------------------------------------------%
%----------                             plotting modules                                 ----------%
%--------------------------------------------------------------------------------------------------%

function plot_chromatic_integration(ci_data, cellnum, clusternum, para)
%
figure('pos',[400 150 1100 750],'color','w');
% first plot chromatic integration curve
subplot(1,2,1)
plot_chromatic_integration_curves(ci_data, para);
% then plot the psths
sp              =   subplot(1,2,2);
plot_chromatic_integration_psths(ci_data, para);
% massive title for both subplots
annotation('textbox', [0.1, 0.95, 0.8, 0], 'string', ['Chromatic integration responses for cell ',...
    num2str(cellnum),', cluster ',num2str(clusternum)],'FontSize',16,'EdgeColor','none',...
    'HorizontalAlignment','center','VerticalAlignment','middle');
% stretch the second plot a bit.
sp.Position     =       [sp.Position(1)-0.05, sp.Position(2)-0.07 sp.Position(3:4)+0.05];
%annotation('textbox', [0.5, 0.06, 0, 0.8], 'string',sprintf('%2g\t\n\n',1:11),'VerticalAlignment','cap');

end

%--------------------------------------------------------------------------------------------------%

function [ci, datplt]= plot_chromatic_integration_curves(ci_data, para)
%
stimpsth        =       ci_data.stimulusPSTH;
pfrpsth         =       ci_data.preframePSTH;

ptch        =   @(x,y,ysem,col)(patch([x,fliplr(x)],[y-ysem,fliplr(y+ysem)],1,...
    'facecolor',col,'edgecolor','none','facealpha',0.3));
pltline     =   @(x,y,col)plot(x, y,'-o','color',col,'markerfacecolor',col,...
    'markeredgecolor',col,'markersize',7,'linewidth',2);
gNuFcol         =       [1 0.0781 0.5742];  % green on, uv off color
gFuNcol         =       [0 0.7461 1];       % green off, uv on color
noaxis          =       false;

contnum         =       size(stimpsth,2)/2;

ci              =       chromatic_integration_curves(stimpsth, pfrpsth, para);
% plotting...plotting...plotting
hold on;
ptch(1 : contnum, ci.grONuvOFF, ci.grONuvOFFsem, gNuFcol);
ptch(1 : contnum, ci.grOFFuvON, ci.grOFFuvONsem, gFuNcol);

respdiff        =       [ci.grONuvOFF+ci.grONuvOFFsem, ci.grONuvOFF-ci.grONuvOFFsem, ...
    ci.grOFFuvON+ci.grOFFuvONsem, ci.grOFFuvON-ci.grOFFuvONsem];

datplt.plot1    =   pltline(1:contnum, ci.grONuvOFF,gNuFcol);
datplt.plot2    =   pltline(1:contnum, ci.grOFFuvON,gFuNcol);
datplt.zeroline =   plot(linspace(0.5,contnum + 0.5,contnum),zeros(1,contnum),'--k','linewidth',0.5);

yAx     =       round(linspace(floor(min(respdiff)/2)*2,ceil(max(respdiff)/2)*2, 4));
if all(yAx == 0),    yAx    =   [0,10];                     end
if all(yAx >  0),    yAx    =   [0,yAx];                    end
if any(yAx <= 0),    yAx    =   unique(sort([0,yAx]));      end


axis([0, 1+contnum, yAx(1)-1, yAx(end)+1]);      axis square;     box off;
yAxneg      =       yAx(yAx<0);    yAxpos   =   yAx(yAx>0);
if not(isempty(yAxneg))
    yAxneg      =   [ceil(yAxneg(1)/5)*5, floor(yAxneg(2:end)/5)*5];
    if all(yAxneg == 0), yAxneg     =   yAx(1);             end
end
yAxpos          =   floor(yAxpos/5)*5;  yAxpos          =       yAxpos(yAxpos~=0);
ylab            =   unique([yAxneg,0,yAxpos]);
xticks(1 : 1 : contnum);        yticks(ylab);
ax = gca;                       ax.TickDir = 'out';

if noaxis
    yrn         =   (yAx(end)+1)- (yAx(1)-1);  %#ok
    yrn         =   ceil((yrn/10)/5)*5;
    if not(ishold),         hold on;            end
    plot([1,1],[ylab(end)/4, yrn + ylab(end)/4],'k','linewidth',1);
    text(1.2, ylab(end)/4 + yrn/2, [num2str(yrn),' Hz'], 'fontsize', 8);
    box off;            ax.YTick  =  [];          ax.YColor  =  'w';
end

xlabel('Stimulus index');
ylabel('Rate (Hz)');
title('Chromatic Integration Curve');

legend([datplt.plot1, datplt.plot2],'Green-On-UV-Off','Green-Off-UV-On','location','southoutside');
legend boxoff;

end

%--------------------------------------------------------------------------------------------------%

function plot_chromatic_integration_psths(ci_data, para)
%
stimpsth        =       ci_data.stimulusPSTH;
pfrpsth         =       ci_data.preframePSTH;
gNuFcol         =       [1 0.0781 0.5742];  % green on, uv off color
gFuNcol         =       [0 0.7461 1];       % green off, uv on color
contnum         =       size(stimpsth,2)/2;
ygap            =       1.15;
xgap            =       0.2;
pfrcol          =       0.75*[1 1 1]; % silver gray color for preframes

% re-arrange input data
cipsth          =       ([pfrpsth(:,1:contnum);stimpsth(:,1:contnum);...
                        pfrpsth(:,1+contnum:2*contnum);stimpsth(:,1+contnum:2*contnum)]);
                    
para.pfrduration    =   para.stimduration;
pfrnumbin           =   para.pfrduration  / para.binlength;
stimnumbin          =   para.stimduration / para.binlength;
timevec             =   [linspace(0,para.pfrduration,pfrnumbin);...
    linspace(para.pfrduration,para.pfrduration+para.stimduration,stimnumbin);...
    linspace(para.pfrduration+para.stimduration,2*para.pfrduration+para.stimduration,pfrnumbin);...
    linspace(2*para.pfrduration+para.stimduration,2*(para.pfrduration+para.stimduration),stimnumbin)];
tIdx                =   reshape(1:length(timevec(:)),[],4)';

ygap                =   ceil((max(cipsth(:))* ygap)/10)*10;
if (ygap == 0),     ygap    =   10;         end
axis([-xgap, (xgap+para.stimduration+para.pfrduration)*2 ,-ygap/4-15, ygap*contnum]);
hold on;        box off;
% add shades
shadefun(xgap);
% draw bars first
for ii = 1 : contnum % start from below
    
    % plotting bars
    barfun(timevec(1,:), cipsth(tIdx(1,:),ii), ((ii-1)* ygap), pfrcol,0.5/2);
    barfun(timevec(2,:), cipsth(tIdx(2,:),ii), ((ii-1)* ygap), gNuFcol,0.5);
    barfun(xgap + timevec(3,:), cipsth(tIdx(3,:),ii), ((ii-1)* ygap),pfrcol,0.5/2);
    barfun(xgap + timevec(4,:), cipsth(tIdx(4,:),ii), ((ii-1)* ygap),gFuNcol,0.5);
    scalebarfun(xgap, ygap, para, ii);
end

stairsfun([timevec(1,:),timevec(2,:)], cipsth(1:size(cipsth,1)/2,:), ygap, contnum, 0, 0.5);
stairsfun([xgap + timevec(3,:),xgap + timevec(4,:)],cipsth(1+size(cipsth,1)/2:size(cipsth,1),:),...
    ygap,contnum, 0, 0.5);
datascalefun(xgap+0.05, ygap, para);
% getting rid of labels
set(gca,'xtick', [],'ytick', [],'xcolor', 'w','ycolor', 'w');
pbaspect([1 1.6 1]);


end

%--------------------------------------------------------------------------------------------------%

function stairsfun(x,y,pltYgap,contnum,drawbaseline,lw)
%
gaps        =       0 : pltYgap : pltYgap * (contnum - 1);
yAx         =       y + repmat(gaps, size(y,1),1);
yAx         =       [gaps ; yAx ; yAx(end,:); gaps];           % this is to close the ending of stairs function
xAx         =       [x(1), x ,  x(end), x(end)];               % this is to make x axis match the ending
stairs( xAx, yAx, 'k', 'linewidth', lw);
if drawbaseline
    if ~ishold,     hold on;        end
    line([xAx(1), xAx(end)],[yAx(1,:) ; yAx(1,:)], 'color', 'k', 'linewidth', lw);
end

end

%--------------------------------------------------------------------------------------------------%

function b = barfun(x,y,g,col,transP,varargin)
%
if iscolumn(y),         y = y';     end
x       =       [x;x];
y       =       [y;y];
xplt    =       x([2:end end]);
yplt    =       y(1:end);
b = patch([xplt, fliplr(xplt)],[yplt + g, fliplr( repmat(g, size(yplt)))], 0 ,...
    'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', transP);
end

%--------------------------------------------------------------------------------------------------%

function shadefun(xgap, shaderegion, col)
%
if nargin <2,       shaderegion = [0.55,0.75,1.55,1.75];        end
if nargin <3,       col = [0.9792 0.9792 0.8203];               end
patch([shaderegion(1), shaderegion(1), shaderegion(2), shaderegion(2)],...
    [get(gca,'YLim'), fliplr(get(gca,'YLim'))], col, 'edgecolor','none');
patch(xgap + [shaderegion(3), shaderegion(3), shaderegion(4),shaderegion(4)],...
    [get(gca,'YLim'), fliplr(get(gca,'YLim'))], col, 'edgecolor','none');

end

%--------------------------------------------------------------------------------------------------%

function scalebarfun(xgap, ygap, para, iter )
%
gr      =       para.contrastgreen;
uv      =       para.contrastuv;
grcol   =       [0 0.8008 0];
uvcol   =       [0.5391 0.1680 0.8828];
tranp   =       0.5;
sc      =       (ygap*0.3)/max(gr);     % scaling ratio for scale bars
tsc     =       (sc*-8)/2;              % text distance from border of each scale bar

barfun([-0.03,0.015], [sc*gr(iter), sc*gr(iter)], ((iter-1)* ygap) + ygap/2, grcol, tranp);
barfun([-0.03,0.015], [sc*uv(iter), sc*uv(iter)], (((iter-1)* ygap) + ygap/2), uvcol, tranp);

text(-0.06,tsc+((((iter-1)* ygap)+ygap/2) + sc*gr(iter)),num2str(gr(iter)),'horiz','right','fontsize',6);
text(-0.06,((((iter-1)* ygap)+ygap/2) + sc*uv(iter))-tsc,num2str(uv(iter)),'horiz','right','fontsize',6);

barfun([xgap+(1.1-0.03),xgap+(1.115)],[-sc*gr(iter),-sc*gr(iter)],((iter-1)* ygap)+ygap/2,grcol,tranp);
barfun([xgap+(1.1-0.03),xgap+(1.115)],[-sc*uv(iter),-sc*uv(iter)],(((iter-1)* ygap)+ygap/2),uvcol,tranp);

text(xgap+(1.1-0.06),((((iter-1)* ygap)+ygap/2) - sc*gr(iter))-tsc,num2str(-gr(iter)),'horiz','right','fontsize',6);
text(xgap+(1.1-0.06),tsc+((((iter-1)* ygap)+ygap/2) - sc*uv(iter)),num2str(-uv(iter)),'horiz','right','fontsize',6);

end

%--------------------------------------------------------------------------------------------------%

function datascalefun(xgap, ygap, para)
% adding scale bar for firing rate
scalloc         =       -ygap/3.5;
tvec            =       (para.stimduration + para.pfrduration) * 2;
scalresp        =       ceil(ygap/50)*10;
plot(xgap + [(0.1/tvec)+tvec tvec tvec],[scalloc scalloc (scalloc + scalresp)],'k-','linewidth',0.5);
text(xgap + tvec + (0.05/tvec),scalloc-5,'100 ms','horiz','center','vert','top','fontsize',7);
text(xgap + tvec - 0.02,(scalloc + scalloc/2)/3,[num2str(scalresp),' spikes'],...
    'horiz','right','vert','middle','fontsize',7);
end

%--------------------------------------------------------------------------------------------------%