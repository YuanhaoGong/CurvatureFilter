function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 12-Jul-2016 22:58:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.axes1,'XTickLabel',[])
set(handles.axes1,'YTickLabel',[])
set(handles.axes2,'XTickLabel',[])
set(handles.axes2,'YTickLabel',[])
set(handles.axes3,'XTickLabel',[])
set(handles.axes3,'YTickLabel',[])

% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
[FileName,PathName] = uigetfile({'*.jpg;*.tif;*.png;','All Image Files';...
 },'Select an image');
global Image;
Image = imread(strcat(PathName,FileName));
global Type;
Type = 2;
set(handles.radiobutton1,'value',0);
set(handles.radiobutton2,'value',0);
set(handles.radiobutton3,'value',1);
set(handles.radiobutton4,'value',0);
global Iteration;
Iteration = 10;
axes(handles.axes1);
imshow(Image);
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Type;
Type = 0;
set(handles.radiobutton1,'value',1);
set(handles.radiobutton2,'value',0);
set(handles.radiobutton3,'value',0);
set(handles.radiobutton4,'value',0);
% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Type;
Type = 1;
set(handles.radiobutton1,'value',0);
set(handles.radiobutton2,'value',1);
set(handles.radiobutton3,'value',0);
set(handles.radiobutton4,'value',0);
% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Type;
Type = 2;
set(handles.radiobutton1,'value',0);
set(handles.radiobutton2,'value',0);
set(handles.radiobutton3,'value',1);
set(handles.radiobutton4,'value',0);
% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Type;
Type = 4;
set(handles.radiobutton1,'value',0);
set(handles.radiobutton2,'value',0);
set(handles.radiobutton3,'value',0);
set(handles.radiobutton4,'value',1);
% Hint: get(hObject,'Value') returns toggle state of radiobutton4



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Iteration;
Iteration = round(str2num(get(hObject,'String')));
if Iteration <= 0 %for some people
    Iteration = - Iteration;
end
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
global Iteration;
Iteration = 10;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Result;
global Image;
global Type;
global Iteration;
global Diff;
Result = CF(Image, Type, Iteration, 1);
Diff = double(Image) - double(Result) + 128;
Diff = uint8(Diff);
Result = uint8(Result);
axes(handles.axes2)
imshow(Result)
axes(handles.axes3)
imshow(Diff)

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Result;
[file,path] = uiputfile('*.png','Save Image');
imwrite(Result,strcat(path,file))

function [result, Energy] = CF(im, FilterType, ItNum, stepsize)
% =========================================================================
% 
%                           Curvature Filter 
% 
% *************************************************************************  
% 
%         @phdthesis{gong:phd, 
%          title={Spectrally regularized surfaces}, 
%          author={Gong, Yuanhao}, 
%          year={2015}, 
%          school={ETH Zurich, Nr. 22616},
%          note={http://dx.doi.org/10.3929/ethz-a-010438292}}
% 
% =========================================================================
% FilterType: 0(Total Variation), 1(Mean Curvature), 2(Gaussian Curvature)
%             4(Bernstein Filter)
if (nargout>1)
    if_record_energy = true;
else
    if_record_energy = false;
end
switch nargin
    case 2
        ItNum = 10; stepsize = 1;
    case 3
        stepsize = 1;
    case 4
    otherwise
        disp('Input parameters are not correct.'), result=0; Energy=0; return;
end

switch FilterType
    case 0
        myfun = @proj_TV; mycurv = @curv_TV;
    case 1
        myfun = @proj_MC; mycurv = @curv_MC;
    case 2
        myfun = @proj_GC; mycurv = @curv_GC;
    case 4
        myfun = @proj_BF; mycurv = @curv_MC;
   otherwise
      disp('Filter Type is not correct.'), return;
end
%pad to even size
orig = single(im); [orig_r, orig_c, orig_z]=size(orig); m=ceil(orig_r/2)*2; n=ceil(orig_c/2)*2;
im = zeros(m,n,orig_z,'single'); im(1:orig_r,1:orig_c,1:orig_z)=orig; 
im(m,:,:) = im(m-1,:,:); im(:,n,:) = im(:,n-1,:); 
%init
result = im; Energy = zeros(ItNum,orig_z);
%% row and column subscript for four types of pixels 
%B = black, W = white, C = circle, T = triangle; r and c indicate row and column
BC_r = uint32(2:2:m-2); BC_c = uint32(2:2:n-2); BT_r = uint32(3:2:m-1); BT_c = uint32(3:2:n-1);
WC_r = uint32(2:2:m-2); WC_c = uint32(3:2:n-1); WT_r = uint32(3:2:m-1); WT_c = uint32(2:2:n-2);
%% the neighbors' index for each type
BC_pre = BC_r-1; BC_nex = BC_r+1; BC_lef = BC_c-1; BC_rig = BC_c+1;
BT_pre = BT_r-1; BT_nex = BT_r+1; BT_lef = BT_c-1; BT_rig = BT_c+1;
WC_pre = WC_r-1; WC_nex = WC_r+1; WC_lef = WC_c-1; WC_rig = WC_c+1;
WT_pre = WT_r-1; WT_nex = WT_r+1; WT_lef = WT_c-1; WT_rig = WT_c+1;
[row,col] = ndgrid(1:size(BT_r,2),1:size(BT_c,2)); row=uint32(row); col=uint32(col);
%% Dual Mesh optimization
for ch = 1:orig_z
    for i = 1:ItNum
        if(if_record_energy) Energy(i, ch) = mycurv(result); %record curvature energy
            if (i>1) && Energy(i,ch) > Energy(i-1,ch) % if the energy start to increase
                break;
            end
        end
        result(:,:,ch) = myfun(result(:,:,ch),BC_r,BC_c,BC_pre,BC_nex,BC_lef,BC_rig,row,col,stepsize);
        result(:,:,ch) = myfun(result(:,:,ch),BT_r,BT_c,BT_pre,BT_nex,BT_lef,BT_rig,row,col,stepsize);
        result(:,:,ch) = myfun(result(:,:,ch),WC_r,WC_c,WC_pre,WC_nex,WC_lef,WC_rig,row,col,stepsize);
        result(:,:,ch) = myfun(result(:,:,ch),WT_r,WT_c,WT_pre,WT_nex,WT_lef,WT_rig,row,col,stepsize);
    end
end
%unpad
result = result(1:orig_r,1:orig_c,1:orig_z);
%% %%%%%%%%%%%%%%%%%%%% three projection operaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = proj_TV(im,BT_r,BT_c,BT_pre,BT_nex,BT_lef,BT_rig,row,col,step)
res = im; BT5 = 5*im(BT_r,BT_c); dist = zeros([size(BT5),8],'single');
tmp1 = im(BT_pre,BT_c) + im(BT_nex,BT_c) - BT5; 
tmp2 = im(BT_r,BT_lef) + im(BT_r,BT_rig) - BT5;
tmp3 = im(BT_pre,BT_c) + im(BT_pre,BT_lef) + im(BT_pre,BT_rig) - BT5; 
tmp4 = im(BT_nex,BT_c) + im(BT_nex,BT_lef) + im(BT_nex,BT_rig) - BT5;
%compute all possible projection distances
dist(:,:,1) = tmp1 + im(BT_pre,BT_lef) + im(BT_r,BT_lef) + im(BT_nex,BT_lef); 
dist(:,:,2) = tmp1 + im(BT_pre,BT_rig) + im(BT_r,BT_rig) + im(BT_nex,BT_rig);
dist(:,:,3) = tmp2 + im(BT_pre,BT_lef) + im(BT_pre,BT_c) + im(BT_pre,BT_rig); 
dist(:,:,4) = tmp2 + im(BT_nex,BT_lef) + im(BT_nex,BT_c) + im(BT_nex,BT_rig); 
dist(:,:,5) = tmp3 + im(BT_r,BT_lef) + im(BT_nex,BT_lef); 
dist(:,:,6) = tmp3 + im(BT_r,BT_rig) + im(BT_nex,BT_rig);
dist(:,:,7) = tmp4 + im(BT_r,BT_lef) + im(BT_pre,BT_lef); 
dist(:,:,8) = tmp4 + im(BT_r,BT_rig) + im(BT_pre,BT_rig);
%find the signed distance with minimal absolute value
tmp = abs(dist); 
[v,ind] = min(tmp,[],3);
%turn sub to index, but faster than sub2ind
dim1 = uint32(size(dist,1));
dim2 = uint32(size(dist,1)*size(dist,2));
index = row + uint32(col-1).*dim1+uint32(ind-1)*dim2;
dm = single(step/5)*dist(index); 
%update current pixels
res(BT_r,BT_c) = res(BT_r,BT_c) + dm;

function res = proj_MC(im,BT_r,BT_c,BT_pre,BT_nex,BT_lef,BT_rig,row,col,step)
res = im; BT8 = 8*im(BT_r,BT_c); dist = zeros([size(BT8),4],'single');
tmp1 = 2.5*(im(BT_pre,BT_c) + im(BT_nex,BT_c)) - BT8;
tmp2 = 2.5*(im(BT_r,BT_lef) + im(BT_r,BT_rig)) - BT8;
%compute all possible projection distances
dist(:,:,1) = tmp1  + 5*im(BT_r,BT_rig) - im(BT_pre,BT_rig) - im(BT_nex,BT_rig);
dist(:,:,2) = tmp1  + 5*im(BT_r,BT_lef) - im(BT_pre,BT_lef) - im(BT_nex,BT_lef);
dist(:,:,3) = tmp2  + 5*im(BT_pre,BT_c) - im(BT_pre,BT_lef) - im(BT_pre,BT_rig);
dist(:,:,4) = tmp2  + 5*im(BT_nex,BT_c) - im(BT_nex,BT_lef) - im(BT_nex,BT_rig);
%find the signed distance with minimal absolute value
tmp = abs(dist); 
[v,ind] = min(tmp,[],3);
%turn sub to index, but faster than sub2ind
dim1 = uint32(size(dist,1));
dim2 = uint32(size(dist,1)*size(dist,2));
index = row + uint32(col-1).*dim1+uint32(ind-1)*dim2;
dm = single(step/8)*dist(index); 
%update current pixels
res(BT_r,BT_c) = res(BT_r,BT_c) + dm;

function res = proj_GC(im,BT_r,BT_c,BT_pre,BT_nex,BT_lef,BT_rig,row,col,step)
res = im; BT2 = 2*im(BT_r,BT_c); BT3 = 3*im(BT_r,BT_c);dist = zeros([size(BT2),8],'single');
%compute all possible projection distances
tmp1 = im(BT_pre,BT_c) - BT3;
tmp2 = im(BT_nex,BT_c) - BT3;
dist(:,:,1) = im(BT_pre,BT_c) + im(BT_nex,BT_c) - BT2; 
dist(:,:,2) = im(BT_r,BT_lef) + im(BT_r,BT_rig) - BT2;
dist(:,:,3) = im(BT_pre,BT_lef) + im(BT_nex,BT_rig) - BT2; 
dist(:,:,4) = im(BT_nex,BT_lef) + im(BT_pre,BT_rig) - BT2;
dist(:,:,5) = tmp1 + im(BT_r,BT_lef) + im(BT_pre,BT_lef); 
dist(:,:,6) = tmp1 + im(BT_r,BT_rig) + im(BT_pre,BT_rig);
dist(:,:,7) = tmp2 + im(BT_r,BT_lef) + im(BT_nex,BT_lef); 
dist(:,:,8) = tmp2 + im(BT_r,BT_rig) + im(BT_nex,BT_rig);
dist(:,:,1:4) = dist(:,:,1:4)*single(1.5); %% scale to the same level
%find the signed distance with minimal absolute value
tmp = abs(dist); 
[v,ind] = min(tmp,[],3); 
%turn sub to index, but faster than sub2ind
dim1 = uint32(size(dist,1));
dim2 = uint32(size(dist,1)*size(dist,2));
index = row + uint32(col-1).*dim1+uint32(ind-1)*dim2;
dm = single(step/3)*dist(index); 
%update current pixels
res(BT_r,BT_c) = res(BT_r,BT_c) + dm;

function res = proj_BF(im,BT_r,BT_c,BT_pre,BT_nex,BT_lef,BT_rig,row,col,step)
res = im; BT2 = 2*im(BT_r,BT_c); BT7 = 7*im(BT_r,BT_c); 
dist = zeros([size(BT2),6],'single');
tmp1 = 3*(im(BT_nex,BT_lef) + im(BT_pre,BT_rig)) - BT7;
tmp2 = 3*(im(BT_pre,BT_lef) + im(BT_nex,BT_rig)) - BT7;
%compute all possible projection distances
dist(:,:,1) = im(BT_pre,BT_c) + im(BT_nex,BT_c) - BT2; 
dist(:,:,2) = im(BT_r,BT_lef) + im(BT_r,BT_rig) - BT2;
dist(:,:,3) = im(BT_pre,BT_c) + im(BT_r,BT_lef) - im(BT_pre,BT_lef) + tmp1; 
dist(:,:,4) = im(BT_pre,BT_c) + im(BT_r,BT_rig) - im(BT_pre,BT_rig) + tmp2;
dist(:,:,5) = im(BT_nex,BT_c) + im(BT_r,BT_lef)- im(BT_nex,BT_lef) + tmp2; 
dist(:,:,6) = im(BT_nex,BT_c) + im(BT_r,BT_rig) - im(BT_nex,BT_rig) +tmp1;
dist(:,:,1:2) = single(3.33333)*dist(:,:,1:2); %% scale to the same level
%find the signed distance with minimal absolute value
tmp = abs(dist); 
[v,ind] = min(tmp,[],3);
%turn sub to index, but faster than sub2ind
dim1 = uint32(size(dist,1));
dim2 = uint32(size(dist,1)*size(dist,2));
index = row + uint32(col-1).*dim1+uint32(ind-1)*dim2;
dm = single(step/10)*dist(index); 
%update current pixels
res(BT_r,BT_c) = res(BT_r,BT_c) + dm;

%% %%%%%%%%%%%%%%%%%%%% curvature energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%these energies are just for observation (so, formula can be changed by user)
function en = curv_TV(im)
[gx,gy]=mygrad(im); 
g = abs(gx) + abs(gy);
en = sum(g(:));
function en = curv_MC(im)
[gx,gy]=mygrad(im);
[gxx,gxy]=mygrad(gx);
[gyx,gyy]=mygrad(gy);
%standard scheme
num = (1+gy.^2).*gxx - gx.*gy.*(gxy+gyx)+ (1+gx.^2).*gyy;
den = (1+gx.^2+gy.^2);
den = sqrt(den).*den*2;
g = num./den;
en = sum(abs(g(:)));
function en = curv_GC(im)
[gx,gy]=mygrad(im);
[gxx,gxy]=mygrad(gx);
[gyx,gyy]=mygrad(gy);
%standard scheme
num = gxx.*gyy-gxy.*gyx;
den = 1+gx.^2+gy.^2;
den = den.*den;
g = num./den; 
en = sum(abs(g(:)));
function [gx, gy]=mygrad(im)
gx=[im(:,2)-im(:,1) (im(:,3:end)-im(:,1:end-2))./2 im(:,end)-im(:,end-1)];
gy=[im(2,:)-im(1,:) ; (im(3:end,:)-im(1:end-2,:))./2 ; im(end,:)-im(end-1,:)];
