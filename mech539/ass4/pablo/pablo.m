function pablo(action);
format short e;
if nargin<1,
    action='initialize';
end;

if strcmp(action,'initialize'),

  global firstcp;
  firstcp = 0;

 back = [0.34 0.67 0.6];
 fore = [0.3 0.3 0.3];

  figNumber=figure( ...
 'Units','normalized', ...
 'Position',[.1 .15 .7 .7],...
 'Name','PABLO', ...
 'NumberTitle','off', ...
 'Color',back, ...
 'Visible','off');
  axes( ...
        'Units','normalized', ...
        'Position',[0.05 0.25 0.7 .70]);
  zoom;

   %====================================
   % Information for all buttons

   xPos=0.77;
   btnWid=0.15;
   btnHt=0.11;
   top=0.35;
   left=0.05;
   right=0.75;
   bottom=0.04;
   spacing=0.04;

   %====================================
   % NBRE OF PANELS 
   btnNumber = 1;
   yPos=0.91-(btnNumber-1)*(btnHt+spacing);
   textStr='  Nb of panels';
   
   % Generic button information
   btnPos1=[xPos yPos btnWid btnHt/3];
   btnPos2=[xPos+btnWid yPos+0.005 btnWid/3 btnHt/3];
   popupHndl=uicontrol( ...
       'Style','text', ...
       'Units','normalized', ... 
       'BackgroundColor',back, ...     
       'ForegroundColor',[0.3 0.3 0.3], ...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ...
       'HorizontalAlignment','left',...
       'FontWeight','bold', ...
       'Position',btnPos1, ...
       'String',textStr);
   callbackStr = 'pablo(''nbpanels'')';
   nbpanelsHndl = uicontrol( ...
	'Style','edit', ...
        'Units','normalized', ...
	'Position',btnPos2, ...
	'Horiz','right', ...
	'Background','white', ...
        'Foreground','black', ...
        'FontUnits','normalized', ...
        'Fontsize',[0.5], ...
	'String','50','Userdata',50, ...
        'callback',callbackStr);

   %====================================
   % The ALFA popup button
   btnNumber=2;
   yPos=1.-(btnNumber-1)*(btnHt+spacing)*.9;
   textStr='            Alpha';
   % Generic button information
   btnPos1=[xPos yPos btnWid btnHt/3];
   btnPos2=[xPos+btnWid yPos+0.005 btnWid/3 btnHt/3];
   alfaHndl=uicontrol( ...
       'Style','text', ...
       'BackgroundColor',back, ...
       'ForegroundColor',[0.3 0.3 0.3], ...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ...
       'HorizontalAlignment','left',...
       'FontWeight','bold', ...
       'Units','normalized', ...
       'Position',btnPos1, ...
      'String',textStr);
   btnPos=[xPos yPos-spacing btnWid btnHt/2];
   callbackStr = 'pablo(''alfa'')';
   alfaHndl = uicontrol( ...
	'Style','edit', ...
        'Units','normalized', ...
	'Position',btnPos2, ...
	'Horiz','right', ...
	'Background','white', ...
        'Foreground','black', ...
        'FontUnits','normalized', ...
        'Fontsize',[.5], ...	
        'String','0','Userdata',0, ...
        'callback',callbackStr);

   %===================================
   % The Viscous radiobutton
   btnNumber = 3;
   yPos=0.96-(btnNumber-1)*(btnHt+spacing)*.7+btnHt/2;
   btnPos=[xPos yPos 1.32*btnWid btnHt/3];
   viscousHndl = uicontrol( ...
       'Style','checkbox', ...
       'Units','normalized', ...
       'Position',btnPos, ...
       'BackgroundColor',back, ...
       'ForegroundColor',fore, ...
       'FontUnits','normalized', ...
       'Value',[0], ...
       'Fontsize',[.5], ...
       'HorizontalAlignment','left',...
       'FontWeight','bold', ...
       'String',' Viscous  ', ...
       'SelectionHighlight','off', ...
       'callback','pablo(''viscousselected'')', ...
       'Interruptible','off');

   textPos=[xPos+0.55*btnWid yPos-btnHt/3+0.002 .77*btnWid btnHt/3];
   callbackStr = 'pablo(''reynolds'')';
   reynoldsHndl = uicontrol( ...
	'Style','edit', ...
	'Units','normalized', ...
	'Position',textPos, ...
	'Horiz','right', ...
	'Background','white', ...
        'Foreground','black', ...
        'FontUnits','normalized', ...
        'Fontsize',[.5], ...
	'String','1000000','Userdata',1000000, ...
	'callback',callbackStr);

    btnPos=[xPos yPos-btnHt/3+0.007 0.52*btnWid btnHt/5];  
    justtext = uicontrol( ...
       'Style','text', ...
       'Units','normalized', ...  
       'BackgroundColor',back, ...
       'ForegroundColor','k', ...
       'Position',btnPos, ...
       'HorizontalAlignment','center',...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ...
       'String','Reynolds :', ...
       'Interruptible','off');  

   %===================================
   % The Ellipse radiobutton
   btnNumber = 3;
   yPos=0.875-(btnNumber-1)*(btnHt+spacing)*.7+btnHt/2-0.04;
   btnPos=[xPos yPos .8*btnWid btnHt/3];
   ellipseHndl = uicontrol( ...
       'Style','checkbox', ...
       'Units','normalized', ...
       'Position',btnPos, ...
       'BackgroundColor',back, ...
       'ForegroundColor',fore, ...
       'FontUnits','normalized', ...
       'Value',[0], ...
       'Fontsize',[.5], ...
       'HorizontalAlignment','left',...
       'FontWeight','bold', ...
       'String',' Ellipse  ', ...
       'SelectionHighlight','off', ...
       'callback','pablo(''ellipseselected'')', ...
       'Interruptible','off');

   textPos=[xPos+0.95*btnWid yPos .38*btnWid btnHt/3];
   callbackStr = 'pablo(''ratio'')';
   ellipseratio = uicontrol( ...
	'Style','edit', ...
	'Units','normalized', ...
	'Position',textPos, ...
	'Horiz','right', ...
	'Background','white', ...
        'Foreground','black', ...
        'FontUnits','normalized', ...
        'Fontsize',[.5], ...
	'String','1','Userdata',1, ...
	'callback',callbackStr);

    btnPos=[xPos+0.142 yPos+0.034 0.4*btnWid btnHt/5];  
    justtext = uicontrol( ...
       'Style','text', ...
       'Units','normalized', ...  
       'BackgroundColor',back, ...
       'ForegroundColor','k', ...
       'Position',btnPos, ...
       'HorizontalAlignment','center',...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ...
       'String','Ratio', ...
       'Interruptible','off');  


   %===================================
   % The NACA editor
   btnNumber = 3;
   yPos=0.8795-(btnNumber-1)*(btnHt+spacing)*.7-0.03;
   btnPos=[xPos yPos .8*btnWid btnHt/3];
   nacaHndl = uicontrol( ...
       'Style','checkbox', ...
       'Units','normalized', ...
       'Position',btnPos, ...
       'BackgroundColor',back, ...
       'ForegroundColor',fore, ...
       'FontUnits','normalized', ...
       'Value',[1], ...
       'Fontsize',[.5], ...
       'HorizontalAlignment','left',...
       'FontWeight','bold', ...
       'String',' NACA  ', ...
       'SelectionHighlight','off', ...
       'callback','pablo(''nacaselected'')', ...
       'Interruptible','off');

   textPos=[xPos+.8*btnWid yPos .16*btnWid btnHt/3];
   callbackStr = 'pablo(''NACA1'')';
   naca1Hndl = uicontrol( ...
	'Style','edit', ...
	'Units','normalized', ...
	'Position',textPos, ...
	'Horiz','right', ...
	'Background','white', ...
        'Foreground','black', ...
        'FontUnits','normalized', ...
        'Fontsize',[.5], ...
	'String','0','Userdata',0, ...
	'callback',callbackStr);

   textPos=[xPos+.8*btnWid+.15*btnWid yPos .16*btnWid btnHt/3];
   callbackStr = 'pablo(''NACA2'')';
   naca2Hndl = uicontrol( ...
	'Style','edit', ...
	'Units','normalized', ...
	'Position',textPos, ...
	'Horiz','right', ...
	'Background','white', ...
	'Foreground','black', ...
        'FontUnits','normalized', ...
        'Fontsize',[.5], ... 
	'String','0','Userdata',0, ...
	'callback',callbackStr);

   textPos=[xPos+.8*btnWid+.3*btnWid yPos .23*btnWid btnHt/3];
   callbackStr = 'pablo(''NACA34'')';
   naca34Hndl = uicontrol( ...
	'Style','edit', ...
	'Units','normalized', ...
	'Position',textPos, ...
	'Horiz','right', ...
	'Background','white', ...
	'Foreground','black', ...
        'FontUnits','normalized', ...
        'Fontsize',[.5], ...
	'String','12','Userdata',12, ...
	'callback',callbackStr);

   %===================================
   % The Airfoil Library button
   btnNumber = 4;
   yPos=0.88-(btnNumber-1)*(btnHt+spacing)*.58-0.02;
   btnPos=[xPos yPos .8*btnWid btnHt/3];
   libHndl = uicontrol( ...
       'Style','checkbox', ...
       'Units','normalized', ...
       'Position',btnPos, ...
       'BackgroundColor',back, ...
       'ForegroundColor',fore, ...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ...
       'HorizontalAlignment','left',...
       'FontWeight','bold', ...
       'SelectionHighlight','off', ...
       'String',' Library  ', ...
       'callback','pablo(''libselected'')', ...
       'Interruptible','off');

   %===================================
   % The Airfoil Library browser
   btnPos=[xPos+ 0.8*btnWid+0.001 yPos .52*btnWid btnHt/3];  
   libpush = uicontrol( ...
       'Style','push', ...
       'Units','normalized', ...
       'Position',btnPos, ...
       'FontUnits','normalized', ...
       'Fontsize',[.4], ...
       'String','Select...', ...
       'callback','pablo(''selectinlib'')', ...
       'Interruptible','off');

   %===================================
   % Display current selected airfoil
    btnPos=[xPos yPos-btnHt/3 1.5*btnWid btnHt/3];  
    currentA = uicontrol( ...
       'Style','text', ...
       'Units','normalized', ...  
       'BackgroundColor',back, ...
       'ForegroundColor','k', ...
       'Position',btnPos, ...
       'HorizontalAlignment','left',...
       'FontUnits','normalized', ...
       'Fontsize',[.4], ...
       'String','Current : no selection', ...
       'Interruptible','off');  

   %===================================
   % Type of singularities
   btnNumber = 5;
   yPos=0.9-(btnNumber-1)*(btnHt+spacing)*.67;
   btnPos=[xPos yPos 1.32*btnWid btnHt/3];
   sHndl = uicontrol( ...
       'Style','checkbox', ...
       'Units','normalized', ...
       'Position',btnPos, ...
       'BackgroundColor',back, ...
       'ForegroundColor','r', ...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ...
       'HorizontalAlignment','left',...
       'FontWeight','bold', ...
       'SelectionHighlight','off', ...
       'String',' Cst Source', ...
       'callback','pablo(''sselected'')', ...
       'Interruptible','off');

   btnNumber = 6;
   yPos=yPos - btnHt/3;
   btnPos=[xPos yPos 1.32*btnWid btnHt/3];
   dHndl = uicontrol( ...
       'Style','checkbox', ...
       'Units','normalized', ...
       'Position',btnPos, ...
       'BackgroundColor',back, ...
       'ForegroundColor','b', ...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ...
       'Value',[1], ...
       'HorizontalAlignment','left',...
       'FontWeight','bold', ...
       'SelectionHighlight','off', ...
       'String',' Cst Doublet', ...
       'callback','pablo(''dselected'')', ...
       'Interruptible','off');

   yPos=yPos - btnHt/3;
   btnPos=[xPos yPos 1.32*btnWid btnHt/3];
   vHndl = uicontrol( ...
       'Style','checkbox', ...
       'Units','normalized', ...
       'Position',btnPos, ...
       'BackgroundColor',back, ...
       'ForegroundColor','g', ...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ...
       'HorizontalAlignment','left',...
       'FontWeight','bold', ...
       'SelectionHighlight','off', ...
       'String',' Lin. Vortex', ...
       'callback','pablo(''vselected'')', ...
       'Interruptible','off');

   %====================================
   % The Hold on  button

   yPos=yPos - 2*btnHt/3;
   btnPos=[xPos yPos 1.32*btnWid btnHt/3];
   holdonHndl = uicontrol( ...
       'Style','checkbox', ...
       'Units','normalized', ...
       'Position',btnPos, ...
       'BackgroundColor',back, ...
       'ForegroundColor','k', ...
       'FontUnits','normalized', ...
       'Fontsize',[.4], ...
       'HorizontalAlignment','left',...
       'SelectionHighlight','off', ...
       'String','  Hold on (Cp)', ...
       'Interruptible','off');

   %====================================
   % The plot bl parameters  button

   yPos=yPos - btnHt/3;
   btnPos=[xPos yPos 1.32*btnWid btnHt/3];
   plotblHndl = uicontrol( ...
       'Style','checkbox', ...
       'Units','normalized', ...
       'Position',btnPos, ...
       'BackgroundColor',back, ...
       'ForegroundColor','k', ...
       'FontUnits','normalized', ...
       'Fontsize',[.4], ...
       'HorizontalAlignment','left',...
       'SelectionHighlight','off', ...
       'String','  Plot BL parameters', ...
       'Interruptible','off');

   %====================================
   % The SHAPE  button

   btnWid=0.19;
   spacing=0.05;
   btnHt2 = 0.6*btnHt;
   shapebtn=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[left bottom btnWid btnHt2], ...
        'FontUnits','normalized', ...
        'Fontsize',[.33], ...
        'String','Shape', ...
        'Interruptible','on', ...
        'Callback','pablo(''shape'')');
 
   %====================================
   % The Cp source  button
   cpbtn=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[left+btnWid+spacing bottom btnWid btnHt2], ...
        'FontUnits','normalized', ...
        'Fontsize',[.33], ...
        'String','Flow solution', ...
        'Interruptible','on', ...
        'Callback','pablo(''cp'')');

   %====================================
   % The INFO button
   infoHndl = uicontrol( ...
	'Style','push', ...
        'Units','normalized', ...
        'Position',[left+2*btnWid+2*spacing bottom btnWid btnHt2], ...
        'FontUnits','normalized', ...
        'Fontsize',[.33], ...
        'String','Info', ...
        'Callback','pablo(''info'')');

   %====================================
   % The CLOSE button
   closeHndl = uicontrol( ...
	'Style','push', ...
        'Units','normalized', ...
        'Position',[0.77 bottom 0.198 btnHt2], ...
        'FontUnits','normalized', ...
        'Fontsize',[.33], ...
        'String','Exit', ...
        'Callback','close all');

   % Now uncover the figure
   set(figNumber,'Visible','on');

   hndlList=[nbpanelsHndl alfaHndl nacaHndl libHndl ... 
             naca1Hndl naca2Hndl naca34Hndl libpush ...
             sHndl dHndl vHndl shapebtn cpbtn infoHndl ...
             closeHndl currentA ellipseHndl ellipseratio holdonHndl ...
             viscousHndl reynoldsHndl plotblHndl];

%    Just to remember : 
%
%    nbpanelsHndl   = hndlList(1);
%    alfaHndl       = hndlList(2);
%    nacaHndl       = hndlList(3);
%    libHndl        = hndlList(4);
%    naca1Hndl      = hndlList(5);
%    naca2Hndl      = hndlList(6);
%    naca34Hndl     = hndlList(7);
%    libpush        = hndlList(8);
%    sHndl          = hndlList(9);
%    dHndl          = hndlList(10);
%    vHdnl          = hndlList(11);
%    shapebtn       = hndlList(12);
%    cpbtn          = hndlList(13);
%    infoHndl       = hndlList(14);
%    closeHndl      = hndlList(15);
%    currentA       = hndlList(16);
%    ellipseHndl    = hndlList(17);
%    ellipseratio   = hndlList(18);
%    holdonHndl     = hndlList(19);
%    viscousHndl    = hndlList(20);
%    reynoldsHndl   = hndlList(21);
%    plotblHndl     = hndlList(22);

   set(figNumber,'Visible','on', 'UserData',hndlList);
    
elseif strcmp(action,'nbpanels'),
    v = get(gco,'UserData');
    s = get(gco,'String');
    vv = eval(s,num2str(v));
    if vv>999, vv=v; end; 
    if vv<0, vv = v; end;
    if rem(vv,2) ~= 0, vv = vv-1; end;
    vv = round(vv);
    set(gco,'Userdata',vv,'String',num2str(vv));

elseif strcmp(action,'alfa'),
    global firstcp;
    firstcp = 0;
    v = get(gco,'UserData');
    s = get(gco,'String');
    vv = eval(s,num2str(v));
    set(gco,'Userdata',vv,'String',num2str(vv));

elseif strcmp(action,'reynolds'),
    mes = 0;
    v = get(gco,'UserData');
    s = get(gco,'String');
    vv = eval(s,num2str(v));
    if vv>40e6, vv=v; mes = 1; end; 
    if vv<1e5, vv = v; mes = 1; end;
    if mes== 1, 
        cla;
        title('Transition model only valid for Reynolds Number between 1e5 and 40e6 !', ...
          'Fontsize',[14],'Color','k');
    end;
    set(gco,'Userdata',vv,'String',num2str(vv));

elseif strcmp(action,'nacaselected'),
    global firstcp;
    firstcp = 0;
    hndlList=get(gcf,'UserData');
    nacaHndl = hndlList(3);
    libHndl = hndlList(4);
    ellipseHndl = hndlList(17);
    if get(nacaHndl,'value')==0
       set(nacaHndl,'value',1);
    else
       set(ellipseHndl,'value',0);
       set(libHndl,'value',0);
    end; 

elseif strcmp(action,'libselected'),
    global firstcp;
    firstcp = 0;
    hndlList=get(gcf,'UserData');
    nacaHndl = hndlList(3);
    libHndl = hndlList(4);
    ellipseHndl = hndlList(17);
    if get(libHndl,'value')==0
       set(libHndl,'value',1);
    else
       set(nacaHndl,'value',0);
       set(ellipseHndl,'value',0);
    end; 

elseif strcmp(action,'ellipseselected'), 
    global firstcp;
    firstcp = 0;
    hndlList=get(gcf,'UserData');
    nacaHndl = hndlList(3);
    libHndl = hndlList(4);
    ellipseHndl = hndlList(17);
    if get(ellipseHndl,'value')==0
       set(ellipseHndl,'value',1);
    else
       set(nacaHndl,'value',0);
       set(libHndl,'value',0);
    end; 

elseif strcmp(action,'sselected'),
    hndlList=get(gcf,'UserData');
    sHndl = hndlList(9);
    dHndl = hndlList(10);
    vHndl = hndlList(11);
    if get(sHndl,'value')==0
       set(sHndl,'value',1);
    else
       set(dHndl,'value',0); 
       set(vHndl,'value',0); 
    end

elseif strcmp(action,'dselected'),
    hndlList=get(gcf,'UserData');
    sHndl = hndlList(9);
    dHndl = hndlList(10);
    vHndl = hndlList(11);
    if get(dHndl,'value')==0
       set(dHndl,'value',1);
    else
       set(sHndl,'value',0); 
       set(vHndl,'value',0); 
    end

elseif strcmp(action,'vselected'),
    hndlList=get(gcf,'UserData');
    sHndl = hndlList(9);
    dHndl = hndlList(10);
    vHndl = hndlList(11);
    if get(vHndl,'value')==0
       set(vHndl,'value',1);
    else
       set(dHndl,'value',0); 
       set(sHndl,'value',0); 
    end

elseif strcmp(action,'selectinlib');
    global firstcp;
    firstcp = 0;
    hndlList=get(gcf,'UserData');
    libHndl = hndlList(4);
    if get(libHndl,'value')==1
      filename = uigetfile('*.DAT','Airfoil Library');
      if filename~=0
         select = filename(1,1:size(filename,2)-4);
         currentA = hndlList(16);
         set(currentA,'String',['Current : ',select]);
      end;
    end;

elseif strcmp(action,'NACA1'),
    global firstcp;
    firstcp = 0;
    hndlList=get(gcf,'UserData');
    holdonHndl     = hndlList(19);
    set(holdonHndl,'value',0);
    v = get(gco,'UserData');
    s = get(gco,'String');
    vv = eval(s,num2str(v));
    if vv<0, vv = v; end;
    if vv>9, vv = v; end;
    vv = round(vv);
    set(gco,'Userdata',vv,'String',num2str(vv));

elseif strcmp(action,'NACA2'),
    global firstcp;
    firstcp = 0;
    hndlList=get(gcf,'UserData');
    v = get(gco,'UserData');
    s = get(gco,'String');
    vv = eval(s,num2str(v));
    if vv<0, vv = v; end;
    if vv>9, vv = v; end;
    vv = round(vv);
    set(gco,'Userdata',vv,'String',num2str(vv));

elseif strcmp(action,'NACA34'),
    global firstcp;
    firstcp = 0;
    hndlList=get(gcf,'UserData');
    holdonHndl     = hndlList(19);
    set(holdonHndl,'value',0);
    v = get(gco,'UserData');
    s = get(gco,'String');
    vv = eval(s,num2str(v));
    if vv<1, vv = v; end;
    if vv>99, vv = 12; end;
    vv = round(vv);
    set(gco,'Userdata',vv,'String',num2str(vv));

elseif strcmp(action,'ratio'),
    global firstcp;
    firstcp = 0;
    hndlList=get(gcf,'UserData');
    v = get(gco,'UserData');
    s = get(gco,'String');
    vv = eval(s,num2str(v));
    if vv<0, vv = v; end;
    if vv>999, vv = v; end;
    set(gco,'Userdata',vv,'String',num2str(vv));

elseif strcmp(action,'shape'),

    back = [0.34 0.67 0.6];
    fore = [0.3 0.3 0.3];
    figNumber=watchon;
    hndlList=get(gcf,'UserData');

    global firstcp;
    firstcp = 0;

    btnPos=[.04 0.11 .95 .11];
    coeffHndl=uicontrol( ...
       'Style','text', ...
       'BackgroundColor',back, ...      
       'ForegroundColor',[0. 0. 0.], ...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ... 
       'Units','normalized', ...
       'Position',btnPos, ...
       'String','');

    btnPos=[.8 0.21 .198 .06];  
    coeffHndl=uicontrol( ...
       'Style','text', ...
       'BackgroundColor',back, ...      
       'ForegroundColor',[0. 0. 0.], ...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ... 
       'Units','normalized', ...
       'Position',btnPos, ...
       'String','');

    nbpanelsHndl   = hndlList(1);
    nacaHndl       = hndlList(3);
    libHndl        = hndlList(4);
    naca1Hndl      = hndlList(5);
    naca2Hndl      = hndlList(6);
    naca34Hndl     = hndlList(7);
    currentA       = hndlList(16);
    ellipseratio   = hndlList(18);
    holdonHndl     = hndlList(19);

    nbp=get(nbpanelsHndl,'UserData');

    if get(nacaHndl,'value') == 1
      eps = get(naca1Hndl,'UserData');  
      p   = get(naca2Hndl,'UserData'); 
      to  = get(naca34Hndl,'UserData'); 
      z  = naca4([eps;p;to],[nbp/2,1]);
      plot (z(:,1),z(:,2),'k');
      grid; axis equal; hold off;
      if to>9
        title(['NACA ',num2str(eps),num2str(p),num2str(to)], ...
          'Fontsize',[14],'Color','k');
      else
        title(['NACA ',num2str(eps),num2str(p),'0',num2str(to)], ...
          'Fontsize',[14],'Color','k');
      end;
    elseif get(libHndl,'value') == 1
      select = get(currentA,'String');
      select = select(1,11:size(select,2));
      if ~strcmp(select,'no selection')
        z = library(nbp/2,select);
        plot (z(:,1),z(:,2),'k');
        grid; axis equal; hold off;
        title([select],'Fontsize',[14],'Color','k');
      else
        cla;
        title('No Airfoil Selected !!', ...
          'Fontsize',[14],'Color','k');
      end;
     else
      axratio = get(ellipseratio,'UserData'); 
      z = ellipse(axratio,nbp/2);
      plot (z(:,1),z(:,2),'k');
      grid; axis equal; hold off;
      title(['Ellipse (Axis Ratio =  ',num2str(axratio),')'], ...
          'Fontsize',[14],'Color','k');
     end;

    watchoff(figNumber);


elseif strcmp(action,'cp'),

    global firstcp;
    firstcp = firstcp + 1;

    back = [0.34 0.67 0.6];
    fore = [0.3 0.3 0.3];
    figNumber=watchon;
    hndlList=get(gcf,'UserData');

    btnPos=[.04 0.11 .95 .11];
    coeffHndl=uicontrol( ...
       'Style','text', ...
       'BackgroundColor',back, ...      
       'ForegroundColor',[0. 0. 0.], ...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ... 
       'Units','normalized', ...
       'Position',btnPos, ...
       'String','');

    btnPos=[.8 0.21 .198 .06];  
    coeffHndl=uicontrol( ...
       'Style','text', ...
       'BackgroundColor',back, ...      
       'ForegroundColor',[0. 0. 0.], ...
       'FontUnits','normalized', ...
       'Fontsize',[.5], ... 
       'Units','normalized', ...
       'Position',btnPos, ...
       'String','');

    nbpanelsHndl   = hndlList(1);
    alfaHndl       = hndlList(2);
    nacaHndl       = hndlList(3);
    libHndl        = hndlList(4);
    naca1Hndl      = hndlList(5);
    naca2Hndl      = hndlList(6);
    naca34Hndl     = hndlList(7);
    sHndl          = hndlList(9);
    dHndl          = hndlList(10);
    currentA       = hndlList(16);
    ellipseratio   = hndlList(18);
    holdonHndl     = hndlList(19);
    viscousHndl    = hndlList(20);    
    reynoldsHndl   = hndlList(21);    
    plotblHndl     = hndlList(22);

    nbp=get(nbpanelsHndl,'UserData');

    alfa = get(alfaHndl,'UserData');
    airfoilselected = 1;

    if firstcp==1
       set(holdonHndl,'value',0);
    end;

    hold off; holdon=0;
    if get(holdonHndl,'value') == 1
      hold on;
      holdon = 1;
      grid;
    end;

    if get(nacaHndl,'value') == 1
      eps = get(naca1Hndl,'UserData'); 
      p   = get(naca2Hndl,'UserData') 
      to  = get(naca34Hndl,'UserData'); 
      z  = naca4([eps;p;to],[nbp/2,1]);

    elseif get(libHndl,'value') == 1
      select = get(currentA,'String');
      select = select(1,11:size(select,2));
      if ~strcmp(select,'no selection')
        z = library(nbp/2,select);
      else        
        airfoilselected = 0;
        cla;
        title('No Airfoil Selected !!', ...
          'Fontsize',[14],'Color','k');
      end;
    else
      axratio = get(ellipseratio,'UserData'); 
      z = ellipse(axratio,nbp/2);
    end;
    
    if airfoilselected ==1
    
       if get(sHndl,'value') == 1 
          clcm = source(z,holdon);
       elseif get(dHndl,'value') == 1 
          clcm = doublet(z,alfa,holdon);
       else
          clcm = vortex(z,alfa,holdon);
       end;
  
       cl = clcm(1);
       cm = clcm(2);
       ue = clcm(3:nbp+3);
       xcp =clcm(nbp+4);
       cpmin = clcm(nbp+5);
       cpmax = clcm(nbp+6);
    
       title('Pressure Coefficient','Fontsize',[14], ...
            'Color','k');   
       grid; hold off;

       if get(viscousHndl,'value') == 1      

          plotbl = 0;
          if get(plotblHndl,'value') == 1, plotbl = 1; end;
               
          % Research of the stagnation point
          isp = stagnation_point(ue);
          ue(isp) = 0;

          % Boundary layer discretization
          nup = isp; nlo = nbp+2-isp;
 
          chordL = 1;
          zup = z(isp:-1:1,:)/chordL;       % Normalized length
          zlo = z(isp:nbp+1,:)/chordL;
   
          Vzero = 1;
          ueup = -ue(isp:-1:1)/Vzero;       % Normalized velocity
          uelo = ue(isp:nbp+1)/Vzero; 

          % Computing the boundary layer's theta and H

          ReL = get(reynoldsHndl,'UserData');

          resup = solvebl(ReL,zup,nup,ueup,plotbl,1);

          figure(figNumber);

          if resup(4)==0
            textStr = ['Upper side : Fully Laminar']; 
          elseif resup(4)==1
            textStr = ['Upper side : Transition at ',num2str(resup(5)),' %']; 
            xt = resup(5)/100;
            yt = findy(z,xt,1)*(cpmax-cpmin);
            h = text(xt,-yt,'T','FontSize',[14],'HorizontalAlignment','center');
          elseif resup(4) ==2
            textStr = ['Upper side : Laminar Separation at ',num2str(resup(5)),' %']; 
            xt = resup(5)/100;
            yt = findy(z,xt,1)*(cpmax-cpmin);
            h = text(xt,-yt,'LS','FontSize',[14],'HorizontalAlignment','center');
          end;
          fileID = fopen('transition.dat','w');
          fprintf(fileID,textStr);
          fprintf(fileID,'\n');

          btnPos=[.05 0.18 .4 .03];
          coeffHndl=uicontrol( ...
            'Style','text', ...
            'BackgroundColor',back, ...
            'ForegroundColor','k', ...
            'HorizontalAlignment','left',...
            'FontUnits','normalized', ...
            'Fontsize',[.6], ...
            'Units','normalized', ...
            'Position',btnPos, ...
            'String',textStr);
              
          if resup(6) ~= 0
            textStr = ['Turbulent Separation at ',num2str(resup(6)),' %']; 
            btnPos=[.455 0.18 .295 .03];
            coeffHndl=uicontrol( ...
              'Style','text', ...
              'BackgroundColor',back, ...
              'ForegroundColor','k', ...
              'HorizontalAlignment','left',...
              'FontUnits','normalized', ...
              'Fontsize',[.6], ...
              'Units','normalized', ...
              'Position',btnPos, ...
              'String',textStr);
            xt = resup(6)/100;
            yt = findy(z,xt,1)*(cpmax-cpmin);
            h = text(xt,-yt,'TS','FontSize',[14],'HorizontalAlignment','center');
          end; 

          reslo = solvebl(ReL,zlo,nlo,uelo,plotbl,2);

          figure(figNumber);
 
          if reslo(4)==0
            textStr = ['Lower side : Fully Laminar']; 
          elseif reslo(4)==1
            textStr = ['Lower side : Transition at ',num2str(reslo(5)),' %']; 
            xt = reslo(5)/100;
            yt = findy(z,xt,2)*(cpmax-cpmin);
            h = text(xt,-yt,'T','FontSize',[14],'HorizontalAlignment','center');
          elseif reslo(4) ==2
            textStr = ['Lower side : Laminar Separation at ',num2str(reslo(5)),' %']; 
            xt = reslo(5)/100;
            yt = findy(z,xt,2)*(cpmax-cpmin);
            h = text(xt,-yt,'LS','FontSize',[14],'HorizontalAlignment','center');
          end;
          fprintf(fileID,textStr);
          fclose(fileID);

          btnPos=[.05 0.14 .4 .03];
          coeffHndl=uicontrol( ...
            'Style','text', ...
            'BackgroundColor',back, ...
            'ForegroundColor','k', ...
            'HorizontalAlignment','left',...
            'FontUnits','normalized', ...
            'Fontsize',[.6], ...
            'Units','normalized', ...
            'Position',btnPos, ...
            'String',textStr);
              
          if reslo(6) ~= 0
            textStr = ['Turbulent Separation at ',num2str(reslo(6)),' %']; 
            btnPos=[.455 0.14 .295 .03];
            coeffHndl=uicontrol( ...
              'Style','text', ...
              'BackgroundColor',back, ...
              'ForegroundColor','k', ...
              'HorizontalAlignment','left',...
              'FontUnits','normalized', ...
              'Fontsize',[.6], ...
              'Units','normalized', ...
              'Position',btnPos, ...
              'String',textStr);
            xt = reslo(6)/100;
            yt = findy(z,xt,2)*(cpmax-cpmin);
            h = text(xt,-yt,'TS','FontSize',[14],'HorizontalAlignment','center');
          end; 

          % Computation of the Drag coefficient

          cd = sy(resup(1),resup(2),resup(3),reslo(1),reslo(2),reslo(3))
          fileID = fopen('cl-cd.dat','a');
          fprintf(fileID,'%f\n',cd);
          fclose(fileID);

       end;
      
       %------ Display the coefficients
       back = [0.34 0.67 0.6];
       fore = [0.3 0.3 0.3];

       btnPos=[.76 0.16 0.24 .11];
       coeffHndl=uicontrol( ...
         'Style','text', ...
         'BackgroundColor',back, ...      
         'ForegroundColor',[0. 0. 0.], ...
         'FontUnits','normalized', ...
         'Fontsize',[.5], ... 
         'Units','normalized', ...
         'Position',btnPos, ...
         'String','');

       textStr = ['Cl   = ',num2str(cl)]; 
       btnPos=[.81 0.23 .198 .03];
       coeffHndl=uicontrol( ...
         'Style','text', ...
         'BackgroundColor',back, ...
         'ForegroundColor','k', ...
         'HorizontalAlignment','left',...
         'FontUnits','normalized', ...
         'Fontsize',[.6], ...
         'Units','normalized', ...
         'Position',btnPos, ...
         'String',textStr);

      textStr = ['Cm = ',num2str(cm)]; 
      btnPos=[.81 0.20 .198 .03];
      coeffHndl=uicontrol( ...
         'Style','text', ...
         'BackgroundColor',back, ...
         'ForegroundColor','k', ...
         'HorizontalAlignment','left',...
         'FontUnits','normalized', ...
         'Fontsize',[.6], ...
         'Units','normalized', ...
         'Position',btnPos, ...
         'String',textStr);
    
      textStr = ['Xcp = ',num2str(xcp)]; 
      btnPos=[.81 0.17 .4 .03];
      coeffHndl=uicontrol( ...
         'Style','text', ...
         'BackgroundColor',back, ...
         'ForegroundColor','k', ...
         'HorizontalAlignment','left',...
         'FontUnits','normalized', ...
         'Fontsize',[.6], ...
         'Units','normalized', ...
         'Position',btnPos, ...
         'String',textStr);
    
      if get(viscousHndl,'value') == 1      
        textStr = ['Cd = ',num2str(cd)]; 
        btnPos=[.81 0.14 .198 .03];
         coeffHndl=uicontrol( ...
          'Style','text', ...
          'BackgroundColor',back, ...
          'ForegroundColor','k', ...
          'HorizontalAlignment','left',...
          'FontUnits','normalized', ...
          'Fontsize',[.6], ...
          'Units','normalized', ...
          'Position',btnPos, ...
          'String',textStr);
      end;       

   end;

   if get(holdonHndl,'value') == 1
      set(holdonHndl,'value',0);
   end;

   watchoff(figNumber);

elseif strcmp(action,'info'),
    ttlStr='';
    hlpStr= ...           
        ['                                                                             '
	 'Panel Methods Theory can be found in:                                        '
         'Katz and Plotkin : Low Speed Aerodynamics, From Wing Theory To Panel Methods.'
         'McGraw-Hill Inc., 1991                                                       '
         '                                                                             ' 
         '* Cst Source  : Constant-Strength Source Distribution (Neumann BC)           '   
         '* Cst Doublet : Constant-Strength Doublet Distribution (Dirichlet BC)        ' 
         '* Lin.Vortex  : Linear-Strength Vortex Distribution (Neumann BC)             '
         '                                                                             '
         '                                                                             '
         'Integral Boundary Layer Theory can be found in :                             '
         'Jack Moran : An introduction to Theoretical and Computational Aerodynamics.  '
         'John Wiley and sons, 1984                                                    '
         '                                                                             '
         '* Laminar boundary layer : Thwaites                                          '
         '* Transition : Michel                                                        '
         '* Turbulent boundary layer : Head                                            ']; 

    info(ttlStr,hlpStr);                
end;    
