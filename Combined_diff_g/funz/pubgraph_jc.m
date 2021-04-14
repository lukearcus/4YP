%publish graph function (makes graphs pretty and publishable)
function pubgraph_jc(fh,FS,LW,CL)
%fh is the handles to the figure containing the graph
%FS is the Font Size 
%LW is the Line Width
%CL is Color of background
if nargin < 4
    CL = [1,1,1];
    if nargin <3 
        LW =1;
        if nargin < 2
            FS = 16;
        end
    end
end 
if nargin < 3
    LW =1;
end
if nargin < 2
    FS = 14;
end

figure(fh)                                         	%pull the figure forwards
axs = findall(fh, 'Type', 'axes');                 	%get the axes on the figure
set(axs, 'FontSize', FS, 'LineWidth', LW,'Color',CL)%make everything on the axis correct
set(findall(fh, 'Type', 'text'), 'FontSize', FS); 	%make all other text correct ,'FontName','Calibri'
set(findall(fh, 'Type', 'line'), 'LineWidth',LW);	%make all other lines correct
if nargin > 3
    set(fh,'Color',CL)                                  %set the figure background color
end

%for each graph on the figure
for i=1:length(axs)
%get the curent axis
    ax=axs(i); 
%     axis(ax,'tight')
%make sure y axis isn't too tight
    v = ylim;                                     	%get current y axis limits
    sc = 0.02*diff(v);                           	%add in 2% each side
    ylim([v(1)-sc , v(2)+sc])                     	%apply it to the graphs
    % v=xlim; sc=0.02*diff(v); xlim([v(1)-sc v(2)+sc])	%make sure x axis isn't too tight
    
%format X and Y axis tickmarks so all numbers are displayed with the same precision
%    YTick = str2num(get(ax,'YTickLabel'));          %get the current Y axis labels
%    XTick = str2num(get(ax,'XTickLabel'));          %get the current X axis labels
    %check what precision the y axis is
%    if ~any(rem(YTick,1)),
%        yfmt='%-5.0f'; 
%    elseif ~any(rem(YTick*10,1)),
%        yfmt='%-5.1f'; 
%    elseif ~any(rem(YTick*100,1)),
%        yfmt='%-5.2f';
%    elseif ~any(rem(YTick*1000,1)),
%        yfmt='%-5.2f';
 %   else
        yfmt='%-5.4f';
%    end
    %check what precision the x axis is
%     if ~any(rem(XTick,1)),
%         xfmt='%-5.0f'; 
%     elseif ~any(rem(XTick*10,1)),
%         xfmt='%-5.1f'; 
%     elseif ~any(rem(XTick*100,1)),
%         xfmt='%-5.2f';
%     elseif ~any(rem(XTick*1000,1)),
%         xfmt='%-5.3f';
%     else
        xfmt='%-5.4f';
%    end
    %set the axis precision
%    set(ax,'YTick',get(ax,'YTick'),'YTickLabel',num2str(YTick,yfmt));%set all the y labels to the same precision
%    set(ax,'XTick',get(ax,'XTick'),'XTickLabel',num2str(XTick,xfmt));%set all the x labels to the same precision
   
end
     tightfig_jc(fh);
end

function hfig = tightfig_jc(hfig)
% tightfig: Alters a figure so that it has the minimum size necessary to
% enclose all axes in the figure without excess space around them.
% 
% Note that tightfig will expand the figure to completely encompass all
% axes if necessary. If any 3D axes are present which have been zoomed,
% tightfig will produce an error, as these cannot easily be dealt with.
% 
% hfig - handle to figure, if not supplied, the current figure will be used
% instead.

    if nargin == 0
        hfig = gcf;
    end

    % There can be an issue with tightfig when the user has been modifying
    % the contnts manually, the code below is an attempt to resolve this,
    % but it has not yet been satisfactorily fixed
%     origwindowstyle = get(hfig, 'WindowStyle');
    set(hfig, 'WindowStyle', 'normal');
    
    % 1 point is 0.3528 mm for future use

    % get all the axes handles note this will also fetch legends and
    % colorbars as well
    hax = findall(hfig, 'type', 'axes');
    
    % get the original axes units, so we can change and reset these again
    % later
    origaxunits = get(hax, 'Units');
    
    % change the axes units to cm
    set(hax, 'Units', 'centimeters');
    
    % get various position parameters of the axes
    if numel(hax) > 1
%         fsize = cell2mat(get(hax, 'FontSize'));
        ti = cell2mat(get(hax,'TightInset'));
        pos = cell2mat(get(hax, 'Position'));
    else
%         fsize = get(hax, 'FontSize');
        ti = get(hax,'TightInset');
        pos = get(hax, 'Position');
    end
    
    % ensure very tiny border so outer box always appears
    ti(ti < 0.1) = 0.15;
    
    % we will check if any 3d axes are zoomed, to do this we will check if
    % they are not being viewed in any of the 2d directions
    views2d = [0,90; 0,0; 90,0];
    
    for i = 1:numel(hax)
        
        set(hax(i), 'LooseInset', ti(i,:));
%         set(hax(i), 'LooseInset', [0,0,0,0]);
        
        % get the current viewing angle of the axes
        [az,el] = view(hax(i));
        
        % determine if the axes are zoomed
        iszoomed = strcmp(get(hax(i), 'CameraViewAngleMode'), 'manual');
        
        % test if we are viewing in 2d mode or a 3d view
        is2d = all(bsxfun(@eq, [az,el], views2d), 2);
               
        if iszoomed && ~any(is2d)
           error('TIGHTFIG:haszoomed3d', 'Cannot make figures containing zoomed 3D axes tight.') 
        end
        
    end
    
    % we will move all the axes down and to the left by the amount
    % necessary to just show the bottom and leftmost axes and labels etc.
    moveleft = min(pos(:,1) - ti(:,1));
    
    movedown = min(pos(:,2) - ti(:,2));
    
    % we will also alter the height and width of the figure to just
    % encompass the topmost and rightmost axes and lables
    figwidth = max(pos(:,1) + pos(:,3) + ti(:,3) - moveleft);
    
    figheight = max(pos(:,2) + pos(:,4) + ti(:,4) - movedown);
    
    % move all the axes
    for i = 1:numel(hax)
        
        set(hax(i), 'Position', [pos(i,1:2) - [moveleft,movedown], pos(i,3:4)]);
        
    end
    
    origfigunits = get(hfig, 'Units');
    
    set(hfig, 'Units', 'centimeters');
    
    % change the size of the figure
    figpos = get(hfig, 'Position');
    
    set(hfig, 'Position', [figpos(1), figpos(2), figwidth, figheight]);
    
    % change the size of the paper
    set(hfig, 'PaperUnits','centimeters');
    set(hfig, 'PaperSize', [figwidth, figheight]);
    set(hfig, 'PaperPositionMode', 'manual');
    set(hfig, 'PaperPosition',[0 0 figwidth figheight]);    
    
    % reset to original units for axes and figure 
    if ~iscell(origaxunits)
        origaxunits = {origaxunits};
    end

    for i = 1:numel(hax)
        set(hax(i), 'Units', origaxunits{i});
    end

    set(hfig, 'Units', origfigunits);
    
%      set(hfig, 'WindowStyle', origwindowstyle);
     
end