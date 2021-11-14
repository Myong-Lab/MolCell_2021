function varargout = maskingGUI(varargin)
% MASKINGGUI MATLAB code for maskingGUI.fig, requires im2mask.m function
%      MASKINGGUI, by itself, creates a new MASKINGGUI or raises the existing
%      singleton*.
%
%      H = MASKINGGUI returns the handle to a new MASKINGGUI or the handle to
%      the existing singleton*.
%
%      MASKINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MASKINGGUI.M with the given input arguments.
%
%      MASKINGGUI('Property','Value',...) creates a new MASKINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before maskingGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to maskingGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help maskingGUI

% Last Modified by GUIDE v2.5 21-Feb-2020 15:47:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @maskingGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @maskingGUI_OutputFcn, ...
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



% --- Executes just before maskingGUI is made visible.
function maskingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to maskingGUI (see VARARGIN)

% Choose default command line output for maskingGUI
handles.output = hObject;

% original image
img = varargin{1};
setappdata(hObject, 'original', img)

% default mask is Mask 6
setappdata(handles.mask_dropdown, 'mask_method', 6)
set(handles.mask_dropdown, 'Value', 6)

% Update handles structure
guidata(hObject, handles);
axes(handles.axes_img)
imshow(imadjust(img))
title('Original')

% UIWAIT makes maskingGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = maskingGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
mask = getappdata(handles.push_chooseMask, 'chosen_mask');
mask_method = getappdata(handles.push_load_mask, 'loaded_mask_method');
varargout{1} = mask;
varargout{2} = mask_method;
close(gcf)


% --- Executes on button press in toggle_orig_mask.
function toggle_orig_mask_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_orig_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_orig_mask

state = get(hObject, 'Value');
mask_method = getappdata(handles.push_load_mask, 'loaded_mask_method');
orig = getappdata(handles.figure1, 'original');
mask = getappdata(handles.push_load_mask, 'loaded_mask');

if state == 0
    axes(handles.axes_img)
    imshow(mask)
    title(['Mask ' num2str(mask_method)])
else
    axes(handles.axes_img)
    imshow(imadjust(orig))
    title('Original')
end


% --- Executes on button press in push_chooseMask.
function push_chooseMask_Callback(hObject, eventdata, handles)
% hObject    handle to push_chooseMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

chosen_mask = getappdata(handles.push_load_mask, 'loaded_mask');
setappdata(handles.push_chooseMask, 'chosen_mask', chosen_mask)
uiresume(handles.figure1)


% --- Executes on selection change in mask_dropdown.
function mask_dropdown_Callback(hObject, eventdata, handles)
% hObject    handle to mask_dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mask_dropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mask_dropdown

mask_method = get(hObject, 'Value');
setappdata(handles.mask_dropdown, 'mask_method', mask_method);


% --- Executes during object creation, after setting all properties.
function mask_dropdown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_load_mask.
function push_load_mask_Callback(hObject, eventdata, handles)
% hObject    handle to push_load_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

img = getappdata(handles.figure1, 'original');
mask_method = getappdata(handles.mask_dropdown, 'mask_method');
setappdata(hObject, 'loaded_mask_method', mask_method);

mask = im2mask(img, mask_method);
setappdata(handles.push_load_mask, 'loaded_mask', mask)
axes(handles.axes_img), imshow(mask)
title(['Mask ' num2str(mask_method)])
