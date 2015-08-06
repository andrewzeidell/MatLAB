function [material, d, index, dlength, dwidth, T, Vd, Vg, I_D] = readSuperDuper(SuperDuperASCIIpath)
%Open the file as fileID
%SuperDuperASCIIpath is the path including the filename
fileID = fopen(SuperDuperASCIIpath)
%Pull the material, dielectric thickness, sample index, length, and
%width of the sample from the .txt file
materialcell = textscan(fileID,'%s %s',1,'HeaderLines',3,'delimiter','\t');
material = materialcell{2}; %eg diF-TES ADT
dielectriccell = textscan(fileID,'%s %s',1,'HeaderLines',4,'delimiter','\t');
dcell = dielectriccell{2};
d = str2double(dcell{1}(1:3))*10^-9; %dielectric thickness in m
indexcell = textscan(fileID,'%s %s %s',1,'HeaderLines',4,'delimiter','\t');
index = indexcell{2};
lengthcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
dlength = lengthcell{2};
widthcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
dwidth = widthcell{2};
tempcell = textscan(fileID,'%s %s %s',1,'HeaderLines',1,'delimiter','\t');
T = str2double(tempcell{2});
Vdcell = textscan(fileID,'%s %s %s',1,'HeaderLines',11,'delimiter','\t');
Vd = str2double(Vdcell{3});
%Pull data from superduper---------------------------------------------
datacell = textscan(fileID,'%f %f %f','HeaderLines',2,'delimiter','\t'); %For Super-Duper change to tab delimited
Vg = datacell{1}; %the Gate-Source voltage is column 1
I_D = datacell{3}; %The Drain Current is column 3
fclose(fileID);
return