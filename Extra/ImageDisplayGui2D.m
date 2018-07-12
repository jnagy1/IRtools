function ImageDisplayGui2D(image, h)
%ImageDisplayGui2D  Auxiliary function for IR Tools
%
% This function implements a GUI with two radio buttons.

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

if nargin == 1
    h = [];
end
%
% Create a figure and axes for the GUI, and display image
%
if isempty(h)
    h = figure('Visible','off');
else
    figure(h), clf
end
imagesc(image), axis off, axis image, colormap(gray), colorbar
 
%
% Create pop-up menu for the color map of the display.
% Default will be hot.
%
%ColorMapPopup = 
uicontrol('Style', 'popup',...
   'String', {'gray','hot','jet','parula','hsv','cool'},...
   'Position', [90 -10 80 50],...
   'Callback', @SetColorMap_Callback);

%
% Create pop-up menu for the scale of the display.
% Default will be no scaling.
%
%ScalePopup = 
uicontrol('Style', 'popup',...
   'String', {'original scale','truncate negative pixels','square root scale','log scale'},...
   'Position', [170 -10 160 50],...
   'Callback', @Scale_Callback);

%
% Create pop-up menu to save image.
% Default will be hot.
%
%SavePopup = 
uicontrol('Style', 'popup',...
   'String', {'save as jpeg','save as eps','help'},...
   'Position', [325 -10 160 50],...
   'Callback', @Save_Callback);
       
          
% 
% Make figure visble after adding all components.
% Note that R2014b allows dot notation for setting and getting, e.g.,
%
%    FigureHandle.Visible = 'on';
%
% but earlier versions do not.  Make sure this works for older releases.
%
set(h, 'Visible', 'on');

%
% Next we have all of the callback functions, which are implemented
% as nested functions (this means all of the above variables are 
% accessible to the nested functions, without having to explicitly
% pass them).
%

c = colorbar;
c.FontSize = 16;

   function SetColorMap_Callback(source,~)
      %
      %  This call back function sets the color map for the image
      %  display.
      %
      val = get(source,'Value');
      maps = get(source,'String');
      %
      %  For R2014b and later, we could implment the above steps as:
      %     val = source.Value;
      %     maps = source.String;
      %
      newmap = maps{val};
      colormap(newmap);
   end

   function Scale_Callback(source,~)
      %
      %  This call back function sets the color map for the image
      %  display.
      %
      val = get(source,'Value');
      switch val
          case 1
              imagesc(image), axis('image', 'off'), colorbar
          case 2
              imagesc(max(image, 0)), axis('image', 'off'), colorbar
          case 3
              imagesc(sqrt(max(image,0))), axis('image', 'off'), colorbar
          case 4
              imagesc(log(max(image,0))), axis('image', 'off'), colorbar
      end
          
   end

function Save_Callback(source,~)
      %
      %  This call back function sets the color map for the image
      %  display.
      %
      val = get(source,'Value');
      if ~isreal(h), h = h.Number; end
      switch val
          case 1
              FileNameString = sprintf('Figure%d.jpg', h);
              MessageString = sprintf('Saving file as ''Figure%d.jpg''',h);
              fprintf(1,[MessageString,'\n']);
              print('-noui',FileNameString,'-djpeg')
          case 2
              FileNameString = sprintf('Figure%d.eps', h);
              MessageString = sprintf('Saving file as ''Figure%d.eps''',h);
              fprintf(1,[MessageString,'\n']);
              print('-noui',FileNameString,'-depsc')
          case 3
              MessageString1 = sprintf('If you want to turn off the colorbar, in the Command Window enter:');
              MessageString2 = sprintf('     colorbar off');
              MessageString3 = sprintf('If you want to save to other formats, or to a specified file name,');
              MessageString4 = sprintf('in the Command Window, enter:');
              MessageString5 = sprintf('     print(''-noui'',MyFileName,FormatFlag)');
              MessageString6 = sprintf('See doc print for more help.');
              fprintf(1,[MessageString1,'\n']);
              fprintf(1,[MessageString2,'\n']);
              fprintf(1,[MessageString3,'\n']);
              fprintf(1,[MessageString4,'\n']);
              fprintf(1,[MessageString5,'\n']);
              fprintf(1,[MessageString6,'\n']);
      end
          
   end

end