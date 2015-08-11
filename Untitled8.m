h = uicontrol('Parent',MuPlot_alpha1,'Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
disp('This will print immediately');
uiwait(MuPlot_alpha1); 
disp('This will print after you click Continue');