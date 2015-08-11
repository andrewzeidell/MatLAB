%SCLC
%KPG 3/14/13

clear all, close all
format long
rootdir = pwd; %Current Directory
% pathdir = 'E:\OE Research';
pathdir = pwd;

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

OriginPath = '\Origin Templates\DeviceData_SCLC.opj'; %The filename of the Origin Template
originObj=actxserver('Origin.ApplicationSI'); %Create COM server with Origin
invoke(originObj,'Execute','doc -mc 1;'); %Make origin session visible
invoke(originObj,'IsModified','false'); %suppress prompt for saving project
invoke(originObj,'Load',strcat(pathdir,OriginPath)); %Load the project (which should be located in
                                                     %the same place as this m file
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
    indexcell = textscan(fileID,'%s %s %s',1,'HeaderLines',8,'delimiter','\t');
    index = indexcell{2};
    lengthcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
    L = str2double(lengthcell{2})*10^-6;
    widthcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
    W = str2double(widthcell{2})*10^-6;
    tempcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
    T = str2double(tempcell{2});
    
    %Overwrite
    W = 50*10^-6;
    L = 120*10^-6;
    t = 50*10^-6;
%     cnsize = length(cellname{a});
%     if cnsize>22
%         Tstring = strcat(cellname{a}(18:20),cellname{a}(24));
%     else
%         Tstring = cellname{a}(18:20);
%     end
%     T = str2double(Tstring(1:3));
    Tstring = tempcell{2};
    
    %Variables for Origin--------------------------------------------------
    WSname = strcat('IV',index,Tstring,num2str(fileID));
    WSTemplate = strcat(pathdir,'\Origin Templates\SCLC.otw');
    
    %Constant variables
    epsr = 3; %relative permittivity for silicon dioxide
    eps0 = 8.854e-12; %Permittivity of free space in F/m
    
    %Pull data from superduper---------------------------------------------
    datacell = textscan(fileID,'%f %f','HeaderLines',4,'delimiter','\t');
    Vd = datacell{1} ;
    I_D = datacell{2};
    
    J = I_D./(t*W);
    logJ = log10(J);
    logV = log10(Vd);
    mu = 10000*8*L^3*abs(J)./(9*epsr*eps0*(Vd.^2));
        
    %Update device worksheets in origin------------------------------------
    invoke(originObj,'CreatePage',2,strcat(WSname{1}),WSTemplate);
    invoke(originObj,'PutWorksheet',WSname{1},[Vd I_D logV J logJ mu]);
    
end

%Message to remind saving--------------------------------------------------
%fileattrib(OriginPath,'-w');
ques = questdlg('Save your Project','MATLAB is Finished','OK','OK');
clear originObj
fclose('all');
cd(rootdir)