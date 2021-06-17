function varargout = SFA_GUI(varargin)
% SFA_GUI MATLAB code for SFA_GUI.fig
%      SFA_GUI, by itself, creates a new SFA_GUI or raises the existing
%      singleton*.
%
%      H = SFA_GUI returns the handle to a new SFA_GUI or the handle to
%      the existing singleton*.
%
%      SFA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SFA_GUI.M with the given input arguments.
%
%      SFA_GUI('Property','Value',...) creates a new SFA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SFA_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SFA_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SFA_GUI

% Last Modified by GUIDE v2.5 24-Sep-2017 18:25:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SFA_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SFA_GUI_OutputFcn, ...
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


% --- Executes just before SFA_GUI is made visible.
function SFA_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SFA_GUI (see VARARGIN)
current_folder=pwd;
addpath(genpath(pwd));
set(handles.Imin,'string','0.3');
set(handles.Imax,'string','1.0');
set(handles.Istep,'string','40');
set(handles.Ip,'string','12.61');
set(handles.fitmin,'string','25');
set(handles.fitmax,'string','60'); %eV
set(handles.filter_thickness,'string','200')
set(handles.fixed_intensity,'string','1')

% Choose default command line output for SFA_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SFA_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SFA_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in filter_name.
function filter_name_Callback(hObject, eventdata, handles)
% hObject    handle to filter_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filter_name contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filter_name


% --- Executes during object creation, after setting all properties.
function filter_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in det_gas.
function det_gas_Callback(hObject, eventdata, handles)
% hObject    handle to det_gas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns det_gas contents as cell array
%        contents{get(hObject,'Value')} returns selected item from det_gas


% --- Executes during object creation, after setting all properties.
function det_gas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to det_gas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Imin_Callback(hObject, eventdata, handles)
% hObject    handle to Imin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Imin as text
%        str2double(get(hObject,'String')) returns contents of Imin as a double


% --- Executes during object creation, after setting all properties.
function Imin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Imin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in run_SFA.
function run_SFA_Callback(hObject, eventdata, handles)
% hObject    handle to run_SFA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fix_int_check_val=get(handles.fix_int_check,'value');
    if fix_int_check_val == 0 % we are NOT fixing the intensity
        det_gas_val=get(handles.det_gas,'value');
        det_gas_list=get(handles.det_gas,'string');
        det_gas=det_gas_list{det_gas_val};
        filter_name_val=get(handles.filter_name,'value');
        filter_name_list=get(handles.filter_name,'string');
        filter_name=filter_name_list{filter_name_val};
        Ip=str2num(get(handles.Ip,'string'));
        Imin=str2num(get(handles.Imin,'string'));
        Imax=str2num(get(handles.Imax,'string'));
        Istep=str2num(get(handles.Istep,'string'));
        fitmin=str2num(get(handles.fitmin,'string'));
        fitmax=str2num(get(handles.fitmax,'string'));
        filter_thickness=str2num(get(handles.filter_thickness,'string'));
        [exp_info, exp_data]=SFA_main_function_v3(Imin,Imax,Istep,Ip,fitmin,fitmax,filter_name,filter_thickness,det_gas,fix_int_check_val);
        %Organization of exp_info and exp_data
        %exp_info: col 1 = origin file name, col 2 = wavelength in m, col 3 =
        %minimizing intensity in units of W/cm^2
        %exp_data: each cell contains an array; col 1 = Integer Harmonic
        %Energies (eV), col 2 = Group delay (fs) col 3 = Group delay error (fs)
        handles.exp_info=exp_info;
        handles.exp_data=exp_data;
        guidata(hObject,handles);
    elseif fix_int_check_val==1 % we are fixing the intensity
        fixed_intensity_val=str2num(get(handles.fixed_intensity,'string'));
        det_gas_val=get(handles.det_gas,'value');
        det_gas_list=get(handles.det_gas,'string');
        det_gas=det_gas_list{det_gas_val};
        filter_name_val=get(handles.filter_name,'value');
        filter_name_list=get(handles.filter_name,'string');
        filter_name=filter_name_list{filter_name_val};
        Ip=str2num(get(handles.Ip,'string'));
        Imin=fixed_intensity_val;%str2num(get(handles.Imin,'string'));
        Imax=fixed_intensity_val;%str2num(get(handles.Imax,'string'));
        Istep=1;%str2num(get(handles.Istep,'string'));
        fitmin=str2num(get(handles.fitmin,'string'));
        fitmax=str2num(get(handles.fitmax,'string'));
        filter_thickness=str2num(get(handles.filter_thickness,'string'));
        [exp_info, exp_data]=SFA_main_function_v3(Imin,Imax,Istep,Ip,fitmin,fitmax,filter_name,filter_thickness,det_gas,fix_int_check_val);
        
    end
    


function Ip_Callback(hObject, eventdata, handles)
% hObject    handle to Ip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ip as text
%        str2double(get(hObject,'String')) returns contents of Ip as a double


% --- Executes during object creation, after setting all properties.
function Ip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Imax_Callback(hObject, eventdata, handles)
% hObject    handle to Imax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Imax as text
%        str2double(get(hObject,'String')) returns contents of Imax as a double


% --- Executes during object creation, after setting all properties.
function Imax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Imax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Istep_Callback(hObject, eventdata, handles)
% hObject    handle to Istep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Istep as text
%        str2double(get(hObject,'String')) returns contents of Istep as a double


% --- Executes during object creation, after setting all properties.
function Istep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Istep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fitmin_Callback(hObject, eventdata, handles)
% hObject    handle to fitmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fitmin as text
%        str2double(get(hObject,'String')) returns contents of fitmin as a double


% --- Executes during object creation, after setting all properties.
function fitmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fitmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fitmax_Callback(hObject, eventdata, handles)
% hObject    handle to fitmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fitmax as text
%        str2double(get(hObject,'String')) returns contents of fitmax as a double


% --- Executes during object creation, after setting all properties.
function fitmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fitmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filter_thickness_Callback(hObject, eventdata, handles)
% hObject    handle to filter_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_thickness as text
%        str2double(get(hObject,'String')) returns contents of filter_thickness as a double


% --- Executes during object creation, after setting all properties.
function filter_thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in close_fig.
function close_fig_Callback(hObject, eventdata, handles)
% hObject    handle to close_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  close_all_but_gui;
  



% --- Executes on button press in export_data.
function export_data_Callback(hObject, eventdata, handles)
% hObject    handle to export_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    exp_info=handles.exp_info;
    exp_data=handles.exp_data;
    det_gas_val=get(handles.det_gas,'value');
    det_gas_list=get(handles.det_gas,'string');
    det_gas=det_gas_list{det_gas_val};
    filter_name_val=get(handles.filter_name,'value');
    filter_name_list=get(handles.filter_name,'string');
    filter_name=filter_name_list{filter_name_val};
    Ip=str2num(get(handles.Ip,'string'));
    Imin=str2num(get(handles.Imin,'string'));
    Imax=str2num(get(handles.Imax,'string'));
    Istep=str2num(get(handles.Istep,'string'));
    fitmin=str2num(get(handles.fitmin,'string'));
    fitmax=str2num(get(handles.fitmax,'string'));
    filter_thickness=str2num(get(handles.filter_thickness,'string'));
    %disp(exp_info);
    %disp(exp_info(1,1));
    %disp(char(exp_info(1,1)));
    %disp(exp_info(1,2));
    %disp(exp_info(1,3));
    %disp(exp_data);
    %disp(transpose(cell2mat(exp_data(1,2))));
    output_data_to_file(exp_info,exp_data,det_gas,filter_name,Ip,Imin,Imax,Istep,fitmin,fitmax,filter_thickness);


% --- Executes on button press in fix_int_check.
function fix_int_check_Callback(hObject, eventdata, handles)
% hObject    handle to fix_int_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fix_int_check



function fixed_intensity_Callback(hObject, eventdata, handles)
% hObject    handle to fixed_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixed_intensity as text
%        str2double(get(hObject,'String')) returns contents of fixed_intensity as a double


% --- Executes during object creation, after setting all properties.
function fixed_intensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixed_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
