 %Pull data from superduper---------------------------------------------
    datacell = textscan(fileID,'%f %f %f','HeaderLines',2,'delimiter','\t'); %For Super-Duper change to tab delimited
    Vg = datacell{1}; %the Gate-Source voltage is column 1
    I_D = datacell{3}; %The Drain Current is column 3