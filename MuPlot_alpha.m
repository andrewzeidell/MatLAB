% MUPLOT - Read in Super Duper ASCII and process Mobility, STS, etc.
%
%
% Description
%     Opens ASCII files output by Super Duper using data collected using the 
%     Agilent Semiconductor analyzer. Allows the user to view data, 
%     modify channel widths, and accept or reject data. Outputs to an
%     Origin template
%
% Inputs ([]s are optional)
%   [string] filePath
%       Relative or absolute path to the .mpl file, with or without '.mpl'
%       extension.
% 
% Outputs
%   [1x1 double] numProfiles
%       An integer indicating the number of profiles read in (up to 60).
%   [1xN struct] mplMeta
%       Where N is equal to numProfiles. A struct array containing
%       information about each profile, corresponding to the columns 
%   [1x1 double] dataSize
%       Number of bins of actual lidar signal.
%   [MxN double] coRaw and crossRaw 
%       Where M is equal to dataSize, and N is equal to numProfiles.
%       Contains the raw lidar returns, one profile in each column.
%   [1x1 double] bgSize
%       Number of bins of background data.
%   [MxN double] coBg, crossBg
%       Where M is equal bgSize, and N is equal to numProfiles.
%       Contains the raw background signals for co-polarized and
%       cross-polarized channels, respectively. Each column
%       corresponds to a column of coRaw and crossRaw, respectively.
%   [M-2 x 1 double] zvec
%       Where M is equal to the number of raw datapoints (usually 1299),
%       thus (M-2) accounts for the 2-bin range & data shift.
% 
%
%
% Author
%   Andrew M. Zeidell <andrewzeidell@gmail.com>
%   Using code written by Katelyn Goetz to process raw data
%   Department of Physics 
%   Wake Forest University
%   Winston Salem, NC, USA
%


function varargout = MuPlot(varargin)
% MUPLOT MATLAB code for MuPlot.fig
%      MUPLOT, by itself, creates a new MUPLOT or raises the existing
%      singleton*.
%
%      H = MUPLOT returns the handle to a new MUPLOT or the handle to
%      the existing singleton*.
%
%      MUPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MUPLOT.M with the given input arguments.
%
%      MUPLOT('Property','Value',...) creates a new MUPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MuPlot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MuPlot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MuPlot

% Last Modified by GUIDE v2.5 28-Jul-2015 09:49:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MuPlot_OpeningFcn, ...
                   'gui_OutputFcn',  @MuPlot_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before MuPlot is made visible.
function MuPlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MuPlot (see VARARGIN)
handles.rootdir = pwd; %Current Directory
handles.cellname =0;
handles.PathName =0;
% handles.pathdir = 'E:\OE Research';
handles.pathdir = pwd;
handles.scrsz = get(0,'ScreenSize');
% Choose default command line output for MuPlot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
    cla;
end
% uiwait();

% UIWAIT makes MuPlot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MuPlot_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in compute_devicedata.
function compute_devicedata_Callback(hObject, eventdata, handles)
% hObject    handle to compute_devicedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
set(gca,'fontsize',14);
cla;
handles.pathdir = get(handles.originpathbox,'String');
guidata(hObject, handles);
if isequal(handles.cellname,0)|| isequal(handles.PathName,0)
    ques = errordlg('No files chosen, please use File->Open to select files','OK');
%     exit %This exits MATLAB if no file is selected.
else
cellname = handles.cellname;
whattype = whos('cellname'); %outputs a structure containing whos data
vartype = whattype.class; %finds the class of 'handles.cellname'
truefalse = strcmp(vartype,'char'); %if only 1 file is selected, the variable type will be a char array
if truefalse==1                     %and truefalse==1, otherwise the class is cell and truefalse==0
    i = 1;
else
    i = length(cellname);
end
%Open Origin---------------------------------------------------------------
OriginPath = '\DeviceData_SC.opj'; %The filename of the Origin Template
originObj=actxserver('Origin.ApplicationSI'); %Create COM server with Origin
invoke(originObj,'Execute','doc -mc 1;'); %Make origin session visible
% invoke(originObj,'Execute','doc -m;'); %Make origin session visible
% invoke(originObj,'IsModified','false'); %suppress prompt for saving project
invoke(originObj,'Load',strcat(handles.pathdir,OriginPath)); %Load the project (which should be located in
                                                     %the same place as this m file
mus = [];
mulins = [];
onoffIs = [];
VTsats = [];
VTlins = [];
SubVTs = [];
traps = [];
nums = [];

%Loop to complete for each device's data file
for a=1:i
    cd(handles.PathName) %makes current directory the one containing the files
    if i==1
        [material, d, index, dlength, dwidth, T, Vd, Vg, I_D] = readSuperDuper([handles.PathName handles.cellname]);
%         fileID = fopen(handles.cellname);
        current_file=handles.cellname;
    else
        [material, d, index, dlength, dwidth, T, Vd, Vg, I_D] = readSuperDuper([handles.PathName handles.cellname{a}]);
%         fileID = fopen(handles.cellname{a});
        current_file=handles.cellname{a};
    end
    
    %Pull the material, dielectric thickness, sample index, length, and
    %width of the sample from the .txt file
%     materialcell = textscan(fileID,'%s %s',1,'HeaderLines',3,'delimiter','\t');
%     material = materialcell{2}; %eg diF-TES ADT
%     dielectriccell = textscan(fileID,'%s %s',1,'HeaderLines',4,'delimiter','\t');
%     dcell = dielectriccell{2};
%     %d = str2double(dcell{1}(1:3))*10^-9; %dielectric thickness in m
%     d = 200*10^-9;
%     indexcell = textscan(fileID,'%s %s %s',1,'HeaderLines',4,'delimiter','\t');
%     index = indexcell{2};
%     %     update the device name in the gui
    set(handles.device_location_display_string,'String',index);
    set(handles.current_file_box,'String',current_file);
%     lengthcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
%     widthcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
    set(handles.customlengthbox,'String',dlength);
    set(handles.customwidthbox,'String',dwidth);
    if get(handles.defaultlengthwidthcheckbox,'Value') == 0
        cd(handles.rootdir)
        menu('Set Custom length and Width before pressing OK','OK');
        cd(handles.PathName);
        %     Lstr = inputdlg(strcat('Enter Length ',index));
        % waits for user to click 'Next', gives opportunity to fill in device data
        % uiwait(handles.nextbutton);

% waitfor(handles.nextbutton,'Value',1);


%     else
    end
    Lstr = get(handles.customlengthbox,'String');
    L = str2double(Lstr)*10^-6;
    %     Wstr = inputdlg(strcat('Enter Width ',index));
    Wstr = get(handles.customwidthbox,'String');
    W = str2double(Wstr)*10^-6;
    
%     L = str2double(lengthcell{2})*10^-6;
   
%     W = str2double(widthcell{2})*10^-6;
%     end
%     tempcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
%     T = str2double(tempcell{2});
%     Vdcell = textscan(fileID,'%s %s %s',1,'HeaderLines',11,'delimiter','\t');
%     Vd = str2double(Vdcell{3});
    
    %Variables for Origin--------------------------------------------------
    WSname = strcat('Mu',index,length,'L',num2str(fileID));
    WSTemplate = strcat(handles.pathdir,'\Device Data WS Template.otw');
    
    %Constant variables: epsr may change later with varying dielectrics -
    %this can be made into a ui or pulled from "superduper"
    epsr = 3.9; %relative permittivity for silicon dioxide
    eps0 = 8.854e-12; %Permittivity of free space in F/m
    C_i = epsr*eps0/d; %Capacitance per unit area in F/m^2
    
    %Pull data from superduper---------------------------------------------
%     datacell = textscan(fileID,'%f %f %f','HeaderLines',2,'delimiter','\t'); %For Super-Duper change to tab delimited
%     Vg = datacell{1}; %the Gate-Source voltage is column 1
%     I_D = datacell{3}; %The Drain Current is column 3
    Idlin = abs(I_D);
    sqI_D = sqrt(abs(I_D)); %Take the square-root
    onoffI = abs(min(I_D)/max(I_D)); %on/off current ratio
    
    %Determine the mostlinear portion of IV--------------------------------
    %This loop steps through the sqI_D vs Vg 10 pts at a time (as in pts
    %1-11, then 2-12, then 3-13, etc.), finds the slope of each portion,
    %and the correlation (R-value) of the data to the fit. The correlations
    %are put into matrix A.
    A = [];
    numpts = 4; %Defines the number of points to try and fit data to
    i2 = length(Vg)-numpts; 
    for b = 1:i2
        g = b+numpts;
        linevals = polyfit(Vg(b:g),sqI_D(b:g),1); %fit a line to 10 points
        fitpts = polyval(linevals,Vg(b:g)); %pull values on the fit
        correlation = corrcoef(sqI_D(b:g),fitpts); %compare "" to the original data
        A = [A;correlation(2)]; %Add the R value (correlation) to the matrix
    end
    
    %This loop finds the indices of the 20 highest correlations in matrix A
    %and uses them to pull the slopes corresponding to this data. The
    %highest slope is used to calculate the mobility.
    SlopeMat = [];
    BMat = [];
    N = [];
    
    for c = 1:20
        bestnum = max(A); %highest correlation value
        n = find(A==bestnum);%index of highest correlation value
        N = [N;n]; %put the index in a matrix to recall later
        p = n+numpts; %p is 10 points greater than n
        linearpts = polyfit(Vg(n:p),sqI_D(n:p),1); %fit a line to the points
        mtry = linearpts(1); %this is the slope
        btry = linearpts(2); %this is the y intercept
        BMat = [BMat;btry]; %put the intercepts in a matrix
        SlopeMat = [SlopeMat;mtry]; %put the slopes in a matrix
        A(n)=0; %set the highest value to 0
    end
    m = max(SlopeMat*-1); %find the maximum slope and use it to calculate mobility
    musat1 = 2*m^2*L/(W*C_i); %mobility in m^2/Vs
    musat = musat1*10000; %mobility in cm^2/Vs
    nind = find(SlopeMat==-m); %find the index of the highest slope
    bint = BMat(nind); %which is also the index of the corresponding intercept
    VTsat = bint/m; %threshold voltage
    n = N(nind); %if you want to plot
    p = n+numpts; %""
    
    %Linear Regime Mobility and Threshold Voltage==========================
    Aa = [];
    i3 = length(Vg)-numpts; 
    for bb = 1:i3
        gg = bb+numpts;
        linevals1 = polyfit(Vg(bb:gg),Idlin(bb:gg),1); %fit a line to # points
        fitpts1 = polyval(linevals1,Vg(bb:gg)); %pull values on the fit
        correlation = corrcoef(Idlin(bb:gg),fitpts1); %compare "" to the original data
        Aa = [Aa;correlation(2)]; %Add the R value (correlation) to the matrix
    end
    
    SlopeMat1 = [];
    BMat1 = [];
    Nn = [];
    
    for cc = 1:20
        bestnum1 = max(Aa); %highest correlation value
        nn = find(Aa==bestnum1);%index of highest correlation value
        Nn = [Nn;nn]; %put the index in a matrix to recall later
        pp = nn+numpts; %p is 10 points greater than n
        linearpts1 = polyfit(Vg(nn:pp),Idlin(nn:pp),1); %fit a line to the points
        mtry1 = linearpts1(1); %this is the slope
        btry1 = linearpts1(2); %this is the y intercept
        BMat1 = [BMat1;btry1]; %put the intercepts in a matrix
        SlopeMat1 = [SlopeMat1;mtry1]; %put the slopes in a matrix
        Aa(nn) = 0; %set the highest value to 0
    end
    m1 = max(SlopeMat1*-1) %find the maximum slope and use it to calculate mobility
    mulin1 = m1*L/(abs(Vd)*W*C_i); %mobility in m^2/Vs
    mulin = mulin1*10000 %mobility in cm^2/Vs
    nind1 = find(SlopeMat1==-m1); %find the index of the highest slope
    bint1 = BMat1(nind1); %which is also the index of the corresponding intercept
    VTlin = bint1/m1; %threshold voltage
    nn = Nn(nind1); %if you want to plot
    pp = nn+numpts; %""
       
    %To find the sub-threshold slope; the method is the same as mobility
    spts = 2;
    logI_D = log10(abs(I_D)); %take the common log of the drain current
    logmMat = []; %empty matrix
    i3 = length(Vg)-spts; %8 points instead of 10
    for j = 1:i3
        k = j+spts;
        linevals2 = polyfit(Vg(j:k),logI_D(j:k),1); %fit a line to 3 points
        logmMat = [logmMat;linevals2(1)];
    end
    maxslope = max(logmMat*-1);
    SubVT = 1/maxslope; %the sub-threshold voltage is the portion of the graph with the highest slope
    aa = find(logmMat==(maxslope*-1)); %to plot
    bb = aa+spts; %to plot 
        
    %Graph to assess program; gives IV curve, regime selected to determine
    %mobility, mu, Vt, subVt, and Ion/Ioff. Highlight through "pause" and
    %press ctrl+R to comment and suppress, ctr+T to uncomment and observe
    %======================================================================
    %musat-----------------------------------------------------------------
%     figure, 
cla;
hold on;
    plot(Vg,sqI_D), plot(Vg(n:p),sqI_D(n:p),'+r');
    axis tight;
    title('Transfer Electrical Characterization, Saturation Regime')
    xlabel('V_{GS} (V)')
    ylabel('\surd(I_D)')
    legend('IV curve','Most Linear Regime','location','NorthEast')
    xt = mean(Vg)+5;
    yt = median(sqI_D);
    textmu = strcat('\mu = ',num2str(musat),' cm^2/Vs');
    text(xt,yt,textmu)
    
    keepques = menu('Is this usable data?','Yes','No','Pick Pts by Hand');
    if keepques==2
%         close all
        continue
    elseif keepques==3
%         figure, 
cla;
        hold on, plot(Vg,sqI_D,'.')
        %temporary 07292015  
        xlim([-20 40]);
        ylim([0 sqI_D(21)])
        title('Choose 5 Points')
        xlabel('V_{GS} (V)')
        ylabel('\surd(I_D)')
        legend('IV curve','location','NorthEast')
        findfigs
        [xchoice,ychoice] = ginput(5);
        linechoice = polyfit(xchoice,ychoice,1);
        m = linechoice(1)*-1;
        b = linechoice(2);
        musat1 = 2*m^2*L/(W*C_i); %mobility in m^2/Vs
        musat = musat1*10000; %mobility in cm^2/Vs
        VTsat = b/m;
    end
    
    %mulin-----------------------------------------------------------------
%     figure,
cla;
hold on, plot(Vg,Idlin), plot(Vg(nn:pp),Idlin(nn:pp),'+r');
    axis tight;
    title('Transfer Electrical Characterization, Linear Regime')
    xlabel('V_{GS} (V)')
    ylabel('I_D')
    legend('IV curve','Most Linear Regime','location','NorthEast')
    xtt = mean(Vg)+5;
    ytt = median(Idlin);
    textmu1 = strcat('\mu = ',num2str(mulin),' cm^2/Vs');
    text(xtt,ytt,textmu1)
    
    keepques = menu('Is this usable data?','Yes','No','Pick Pts by Hand');
    if keepques==2
%         close all
        continue
    elseif keepques==3
        cla;
%         figure,
        hold on, plot(Vg,Idlin,'.')
        %temporary 07292015  
        xlim([-20 40]);
        ylim([0 Idlin(21)])
        title('Choose 5 Points')
        xlabel('V_{GS} (V)')
        ylabel('\surd(I_D)')
        legend('IV curve','location','NorthEast')
        findfigs
        [xchoice1,ychoice1] = ginput(5);
        linechoice1 = polyfit(xchoice1,ychoice1,1);
        m1 = linechoice1(1)*-1;
        b1 = linechoice1(2);
        mulin1 = m1*L/(abs(Vd)*W*C_i); %mobility in m^2/Vs
        mulin = mulin1*10000; %mobility in cm^2/Vs
        VTlin = b1/m1;
    end
    
    %Subthreshold Slope----------------------------------------------------   
%     figure, 
cla;
hold on, plot(Vg,logI_D), plot(Vg(aa:bb),logI_D(aa:bb),'+r');
    axis tight;
    title('Sub Threshold Swing')
    xlabel('V_{GS} (V)')
    ylabel('log(I_D)')
    keepques1 = menu('Is this usable data?','Yes','No','Pick Pts by Hand');
    if keepques1==2
%         close all
        continue
    elseif keepques1==3
%         figure,
cla;
hold on, plot(Vg,logI_D,'.')
%temporary 07292015  
        xlim([-20 40]);
        ylim([min(logI_D)  logI_D(21)])
        title('Choose 5 Points')
        xlabel('V_{GS} (V)')
        ylabel('log(I_D)')
        legend('Subthreshold Swing','location','NorthEast')
        findfigs
        [xchoice1,ychoice1] = ginput(5);
        linechoice = polyfit(xchoice1,ychoice1,1);
        aa = linechoice(1)*-1;
        bb = linechoice(2);
        SubVT = 1/aa;
    end
%     close all
    
    %Update device worksheets in origin------------------------------------
    invoke(originObj,'CreatePage',2,strcat(WSname{1}),WSTemplate);
    invoke(originObj,'PutWorksheet',WSname{1},[Vg I_D]);
    %Temp, andrew 07292015
invoke(originObj, 'Execute', ['plotxy [' WSname{1} ']Sheet1!(1,3) plot:=201 ogl:=[<new template:=XY>]']);
invoke(originObj, 'Execute', 'layadd type:=rightY');
invoke(originObj, 'Execute', 'layer.y.type = 2');
invoke(originObj, 'Execute', ['plotxy [' WSname{1} ']Sheet1!(1,4) plot:=201 ogl:=2']);
invoke(originObj, 'Execute', 'label -yr log(Id)');
invoke(originObj, 'Execute', 'yr.label.rotate=270');
invoke(originObj, 'Execute', ['page.longname$ =' WSname{1} 'graph']);
    mus = [mus;musat];
    mulins = [mulins;mulin];
    VTsats = [VTsats;VTsat];
    VTlins = [VTlins;VTlin];
    SubVTs = [SubVTs;SubVT];
    N_i = (1.602e-19*SubVT/(1.38e-23*294*log(10))-1)*C_i/(1.602e-15);
    traps = [traps;N_i];
    nums = [nums;fileID];
    onoffIs = [onoffIs;onoffI];
%     cla;
end

data4origin = [nums,mus,mulins,onoffIs,VTsats,VTlins,SubVTs,traps];

WSTemplatePath = handles.pathdir;
WSTemplate1 = strcat(WSTemplatePath,'\SingleCrystalData.otw');
invoke(originObj,'CreatePage',2,'SCData',WSTemplate1);
invoke(originObj,'PutWorksheet','SCData',data4origin);

%Message to remind saving------------------------------------------------------
%fileattrib(OriginPath,'-w');
ques = warndlg('Save your Project','MATLAB is Finished');
clear originObj
fclose('all');
cd(handles.rootdir)
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cellname,PathName] = uigetfile('*.*','Open the existing ASCII file',...
    'MultiSelect','on');

if isequal(cellname,0)||isequal(PathName,0)
%     exit %This exits MATLAB if no file is selected.
 ques = errordlg('No files chosen, please use File->Open to select files','OK');
end


handles.cellname = cellname;
handles.PathName = PathName;
guidata(hObject, handles);

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function originpathbox_Callback(hObject, eventdata, handles)
% hObject    handle to originpathbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of originpathbox as text
%        str2double(get(hObject,'String')) returns contents of originpathbox as a double


% --- Executes during object creation, after setting all properties.
function originpathbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to originpathbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_rootdir.
function browse_rootdir_Callback(hObject, eventdata, handles)
% hObject    handle to browse_rootdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in choose_originpath.
function choose_originpath_Callback(hObject, eventdata, handles)
% hObject    handle to choose_originpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[InPathName] = uigetdir(get(handles.originpathbox,'String'));
if InPathName ~= 0
    %     InPathName = strcat(InPathName,'/');
    set(handles.originpathbox,'String',[InPathName]);
    handles.pathdir = InPathName;
    guidata(hObject, handles);
end


% --- Executes on button press in defaultlengthwidthcheckbox.
function defaultlengthwidthcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to defaultlengthwidthcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of defaultlengthwidthcheckbox



function customlengthbox_Callback(hObject, eventdata, handles)
% hObject    handle to customlengthbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of customlengthbox as text
%        str2double(get(hObject,'String')) returns contents of customlengthbox as a double


% --- Executes during object creation, after setting all properties.
function customlengthbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to customlengthbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function customwidthbox_Callback(hObject, eventdata, handles)
% hObject    handle to customwidthbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of customwidthbox as text
%        str2double(get(hObject,'String')) returns contents of customwidthbox as a double


% --- Executes during object creation, after setting all properties.
function customwidthbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to customwidthbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nextbutton.
function nextbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nextbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pause off;
