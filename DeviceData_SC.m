%Device Data Calculation for Single Crystal Devices (no stats)

%This program calculates the mobility of an OFET based on the saturation
%regime model: mu = 2*(d(sqrt(I_Dsat))/d(V_G))^2*L/(W*C_i). It also
%calculates subthreshold slope, trap density at the
%dielectric/semiconductor interface, on/off current ratio, and threshold
%voltage

%Katelyn Goetz 
%12/02/11
%2/8/2012
%9/5/2012

clear all, close all
format long
rootdir = pwd; %Current Directory
% pathdir = 'E:\OE Research';
pathdir = 'M:\Matlab\Origin Templates\';
scrsz = get(0,'ScreenSize');

%Variable List-------------------------------------------------------------
%mu  = charge carrier mobility
%C_i = Capacitance per unit area
%L   = crystal length
%W   = crystal width
%I_D = drain current in the saturation regime
%Vg  = gate voltage

%Opens files and checks existance
[cellname,PathName] = uigetfile('*.*','Open the existing ASCII file',...
    'MultiSelect','on');
if isequal(cellname,0)||isequal(PathName,0)
    exit %This exits MATLAB if no file is selected.
end

whattype = whos('cellname'); %outputs a structure containing whos data
vartype = whattype.class; %finds the class of 'cellname'
truefalse = strcmp(vartype,'char'); %if only 1 file is selected, the variable type will be a char array
if truefalse==1                     %and truefalse==1, otherwise the class is cell and truefalse==0
    i = 1;
else
    i = length(cellname);
end

%Open Origin---------------------------------------------------------------
OriginPath = '\Origin Templates\DeviceData_SC.opj'; %The filename of the Origin Template
originObj=actxserver('Origin.ApplicationSI'); %Create COM server with Origin
invoke(originObj,'Execute','doc -mc 1;'); %Make origin session visible
invoke(originObj,'IsModified','false'); %suppress prompt for saving project
invoke(originObj,'Load',strcat(pathdir,OriginPath)); %Load the project (which should be located in
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
    cd(PathName) %makes current directory the one containing the files
    if i==1
        fileID = fopen(cellname);
    else
        fileID = fopen(cellname{a});
    end
    
    %Pull the material, dielectric thickness, sample index, length, and
    %width of the sample from the .txt file
    materialcell = textscan(fileID,'%s %s',1,'HeaderLines',3,'delimiter','\t');
    material = materialcell{2}; %eg diF-TES ADT
    dielectriccell = textscan(fileID,'%s %s',1,'HeaderLines',4,'delimiter','\t');
    dcell = dielectriccell{2};
    %d = str2double(dcell{1}(1:3))*10^-9; %dielectric thickness in m
    d = 200*10^-9;
    indexcell = textscan(fileID,'%s %s %s',1,'HeaderLines',4,'delimiter','\t');
    index = indexcell{2};
    lengthcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
    %L = str2double(lengthcell{2})*10^-6;
    Lstr = inputdlg(strcat('Enter Length ',index));
    L = str2double(Lstr)*10^-6;
    widthcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
    %W = str2double(widthcell{2})*10^-6;
    Wstr = inputdlg(strcat('Enter Width ',index));
    W = str2double(Wstr)*10^-6;
    tempcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
    T = str2double(tempcell{2});
    Vdcell = textscan(fileID,'%s %s %s',1,'HeaderLines',11,'delimiter','\t');
    Vd = str2double(Vdcell{3});
    
    %Variables for Origin--------------------------------------------------
    WSname = strcat('Mu',index,lengthcell{2},'L',num2str(fileID));
    WSTemplate = strcat(pathdir,'\Origin Templates\Device Data WS Template.otw');
    
    %Constant variables: epsr may change later with varying dielectrics -
    %this can be made into a ui or pulled from "superduper"
    epsr = 3.9; %relative permittivity for silicon dioxide
    eps0 = 8.854e-12; %Permittivity of free space in F/m
    C_i = epsr*eps0/d; %Capacitance per unit area in F/m^2
    
    %Pull data from superduper---------------------------------------------
    datacell = textscan(fileID,'%f %f %f','HeaderLines',2,'delimiter','\t'); %For Super-Duper change to tab delimited
    Vg = datacell{1}; %the Gate-Source voltage is column 1
    I_D = datacell{3}; %The Drain Current is column 3
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
    figure, hold on, plot(Vg,sqI_D), plot(Vg(n:p),sqI_D(n:p),'+r')
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
        close all
        continue
    elseif keepques==3
        figure, hold on, plot(Vg,sqI_D,'.')
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
    figure, hold on, plot(Vg,Idlin), plot(Vg(nn:pp),Idlin(nn:pp),'+r')
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
        close all
        continue
    elseif keepques==3
        figure, hold on, plot(Vg,Idlin,'.')
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
    figure, hold on, plot(Vg,logI_D), plot(Vg(aa:bb),logI_D(aa:bb),'+r')
    title('Sub Threshold Swing')
    xlabel('V_{GS} (V)')
    ylabel('log(I_D)')
    keepques1 = menu('Is this usable data?','Yes','No','Pick Pts by Hand');
    if keepques1==2
        close all
        continue
    elseif keepques1==3
        figure, hold on, plot(Vg,logI_D,'.')
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
    close all
    
    %Update device worksheets in origin------------------------------------
    invoke(originObj,'CreatePage',2,strcat(WSname{1}),WSTemplate);
    invoke(originObj,'PutWorksheet',WSname{1},[Vg I_D]);
    
    mus = [mus;musat];
    mulins = [mulins;mulin];
    VTsats = [VTsats;VTsat];
    VTlins = [VTlins;VTlin];
    SubVTs = [SubVTs;SubVT];
    N_i = (1.602e-19*SubVT/(1.38e-23*294*log(10))-1)*C_i/(1.602e-15);
    traps = [traps;N_i];
    nums = [nums;fileID];
    onoffIs = [onoffIs;onoffI];
end

data4origin = [nums,mus,mulins,onoffIs,VTsats,VTlins,SubVTs,traps];

WSTemplatePath = 'E:\OE Research\Origin Templates';
WSTemplate1 = strcat(WSTemplatePath,'\SingleCrystalData.otw');
invoke(originObj,'CreatePage',2,'SCData',WSTemplate1);
invoke(originObj,'PutWorksheet','SCData',data4origin);

%Message to remind saving------------------------------------------------------
%fileattrib(OriginPath,'-w');
ques = questdlg('Save your Project','MATLAB is Finished','OK','OK');
clear originObj
fclose('all');
cd(rootdir)