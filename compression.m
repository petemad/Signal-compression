function varargout = finaltask2(varargin)
% FINALTASK2 MATLAB code for finaltask2.fig
%      FINALTASK2, by itself, creates a new FINALTASK2 or raises the existing
%      singleton*.
%
%      H = FINALTASK2 returns the handle to a new FINALTASK2 or the handle to
%      the existing singleton*.
%
%      FINALTASK2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINALTASK2.M with the given input arguments.
%
%      FINALTASK2('Property','Value',...) creates a new FINALTASK2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before finaltask2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to finaltask2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help finaltask2

% Last Modified by GUIDE v2.5 14-Mar-2019 00:17:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @finaltask2_OpeningFcn, ...
                   'gui_OutputFcn',  @finaltask2_OutputFcn, ...
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


% --- Executes just before finaltask2 is made visible.
function finaltask2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to finaltask2 (see VARARGIN)

% Choose default command line output for finaltask2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes finaltask2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = finaltask2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in compress.
function compress_Callback(hObject, eventdata, handles)
global BFflag
global Lflag 
global signal
switch BFflag
case 0 %DCT
    switch Lflag
        case 0 %Lossy
            dct_Sig=dct(signal);
            Pos_th=find(abs(dct_Sig)<0.5);
            Sig_Afterthreshold=dct_Sig;
            Sig_Afterthreshold(Pos_th)=0;
            
            Q=8;% 1 byte (value representation)
            Max_SigAfter=max(Sig_Afterthreshold);
            Min_SigAfter=min(Sig_Afterthreshold);
            Sig_Quant=round((-1+2^Q)*(Sig_Afterthreshold-Min_SigAfter)/(Max_SigAfter-Min_SigAfter));
            
            Sig_RLE=rle(Sig_Quant);
            save('Sig_RLE');
            
            sigRLE= whos('Sig_RLE');
            sig = whos('signal');
            sig_size = sig.bytes;
            sigRLE_size = sigRLE.bytes;
            CR = sig_size/sigRLE_size;
            set(handles.text2, 'String', CR);
            
        case 1 %LossLess
            dct_Sig=dct(signal);
            
            Sig_RLE=rle(dct_Sig);
            save('Sig_RLE');
            sigRLE= whos('Sig_RLE');
            sig = whos('signal');
            sig_size = sig.bytes;
            sigRLE_size = sigRLE.bytes;

            CR = sig_size/sigRLE_size
            set(handles.text2, 'String', CR); 
    end
case 1 %DFT
    switch Lflag
        case 0 %Lossy
            Sig_fft = fft(signal);
            len=1:1:length(Sig_fft);
            for i = 1:length(len) 
                if abs(Sig_fft(i)) < 300
                    Sig_fft(i) = 0;
                end
            end

            Sig_RLE=rle(Sig_fft);
            save('Sig_RLE');
            
            sigRLE= whos('Sig_RLE');
            sig = whos('signal');
            sig_size = sig.bytes;
            sigRLE_size = sigRLE.bytes;
            CR = sig_size/sigRLE_size;
            set(handles.text2, 'String', CR);

        case 1 %lossless
            Sig_fft = fft(signal);
            
            Sig_RLE=rle(Sig_fft);
            save('Sig_RLE');
            sigRLE= whos('Sig_RLE');
            sig = whos('signal');
            sig_size = sig.bytes;
            sigRLE_size = sigRLE.bytes;
            CR = sig_size/sigRLE_size;
            set(handles.text2, 'String', CR);
    end
case 2 %W1
    switch Lflag
        case 0 %Lossy
            wht_Sig=fwht(signal);
            Pos_th=find(abs(wht_Sig)<0.1);
            Sig_Afterthreshold=wht_Sig;
            Sig_Afterthreshold(Pos_th)=0;
            
            Q=8;% 1 byte (value representation)
            Max_SigAfter=max(Sig_Afterthreshold);
            Min_SigAfter=min(Sig_Afterthreshold);
            Sig_Quant=round((-1+2^Q)*(Sig_Afterthreshold-Min_SigAfter)/(Max_SigAfter-Min_SigAfter));
            
            Sig_RLE=rle(Sig_Quant);
            save('Sig_RLE');
            
            sigRLE= whos('Sig_RLE');
            sig = whos('signal');
            sig_size = sig.bytes;
            sigRLE_size = sigRLE.bytes;
            CR = sig_size/sigRLE_size;
            set(handles.text2, 'String', CR); 
        case 1 %Lossless
            wht_Sig=fwht(signal);
            
            Sig_RLE=rle(wht_Sig);
            Sig_RLE=single(Sig_RLE);
            save('Sig_RLE');
            
            sigRLE= whos('Sig_RLE');
            sig = whos('signal');
            sig_size = sig.bytes;
            sigRLE_size = sigRLE.bytes;
            CR = sig_size/sigRLE_size;
            set(handles.text2, 'String', CR);
    end
case 3 %W2
    switch Lflag
        case 0 %Lossy
            y=1:1:length(signal);

            DEC = mdwtdec('c',signal,5,'sym4');

            [XC,~] = mswcmp('cmp',DEC,'bal_sn');
            len=1:1:length(XC);
            for i = 1:length(len) 
                if abs(XC(i)) < 0.1
                    XC(i) = 0;
                end
            end
            mu = 0;
            sigma = 1;
            pd = makedist('Normal','mu',mu,'sigma',sigma);
            p = pdf (pd,XC);
                p1=p/sum(p);
            XCu = unique(XC);
            len=1:1:length(XCu);
            prob = 1:length(len);
            counter = 0;
            for i = 1:length(len)
                for j = 1:length(y)

                        if (XCu(i)==XC(j))
                          counter = counter+1;
                        end         
                end
                prob(i)=counter/length(y);
                counter = 0;
            end
            [dict,~] = huffmandict(XCu,prob);
            enc = huffmanenco(XC,dict);
            enc = uint8(enc)
            save('enc')
            
            sigEnc= whos('enc');
            sig = whos('signal');
            sig_size = sig.bytes;
            sigEnc_size = sigEnc.bytes;
            CR = sig_size/sigEnc_size;
            set(handles.text2, 'String', CR);
            


            legend('signal ','Diff');
        case 1 %Lossless
            y=1:1:length(signal);

            DEC = mdwtdec('c',signal,5,'sym4');

            [XC,DECCMP,THRESH] = mswcmp('cmp',DEC,'bal_sn');
            mu = 0;
            sigma = 1;
            pd = makedist('Normal','mu',mu,'sigma',sigma);
            p = pdf (pd,XC);
                p1=p/sum(p);
            XCu = unique(XC);
            len=1:1:length(XCu);
            prob = 1:length(len);
            counter = 0;
            for i = 1:length(len)
                for j = 1:length(y)

                        if (XCu(i)==XC(j))
                          counter = counter+1;
                        end         
                end
                prob(i)=counter/length(y);
                counter = 0;
            end
            [dict,~] = huffmandict(XCu,prob);
            enc = huffmanenco(XC,dict);
            enc = uint8(enc);

            save('enc');
            cmpSize=whos('enc');

            cmpSize=cmpSize.bytes;
            OrigSize=whos('signal')
            OrigSize=OrigSize.bytes;
            CR=(OrigSize/cmpSize);
            set(handles.text2, 'String', CR);
            axes(handles.axes2);
            plot ([XC;signal-XC]');  


            legend('XC','Diff');  
    end
end
    


% --- Executes on button press in uncompress.
function uncompress_Callback(hObject, eventdata, handles)
global BFflag
global Lflag 
global signal
switch BFflag
case 0 %DCT
    switch Lflag
        case 0 %lossy
            dct_Sig=dct(signal);
            Pos_th=find(abs(dct_Sig)<0.5);
            Sig_Afterthreshold=dct_Sig;
            Sig_Afterthreshold(Pos_th)=0;
            
            Q=8;% 1 byte (value representation)
            Max_SigAfter=max(Sig_Afterthreshold);
            Min_SigAfter=min(Sig_Afterthreshold);
            Sig_Quant=round((-1+2^Q)*(Sig_Afterthreshold-Min_SigAfter)/(Max_SigAfter-Min_SigAfter));
            
            Sig_RLE=rle(Sig_Quant);
            
            Sig_Irle=irle(Sig_RLE);% Inverse Run Length Codising 
            Sig_IQuant=((Max_SigAfter-Min_SigAfter)/(-1+2^Q))*Sig_Irle+Min_SigAfter;% Inverse Quantization
            Rec_Sig=idct(Sig_IQuant);%Décompressed Signal
            save('Rec_Sig');
            
            axes(handles.axes2); 
            plot([signal-Rec_Sig; Rec_Sig]');
        case 1 %lossless
            dct_Sig=dct(signal);
            Sig_RLE=rle(dct_Sig);
            
            Sig_Irle=irle(Sig_RLE);% Inverse Run Length Codising 
            Rec_Sig=idct(Sig_Irle);%Décompressed Signal
            save('Rec_Sig');
            
            axes(handles.axes2); 
            plot([signal-Rec_Sig; Rec_Sig]');
    end
case 1 %DFT
    switch Lflag 
        case 0 %lossy
            Sig_fft = fft(signal);
            len=1:1:length(Sig_fft);
            for i = 1:length(len) 
                if abs(Sig_fft(i)) < 300
                    Sig_fft(i) = 0;
                end
            end
            
            Sig_RLE=rle(Sig_fft);
            Sig_Irle=irle(Sig_RLE);   
            Rec_Sig=ifft(Sig_Irle);
            save('Rec_Sig');
            
            axes(handles.axes2); 
            plot(Rec_Sig);
        case 1 %lossless
            Sig_fft = fft(signal);
            
            Sig_RLE=rle(Sig_fft);
            Sig_Irle=irle(Sig_RLE);   
            Rec_Sig=ifft(Sig_Irle);
            save('Rec_Sig');
            
            axes(handles.axes2); 
            plot([signal-Rec_Sig; Rec_Sig]');
    end
case 2 %WHT
    switch Lflag
        case 0 %Lossy
            wht_Sig=fwht(signal);
            Pos_th=find(abs(wht_Sig)<0.1);
            Sig_Afterthreshold=wht_Sig;
            Sig_Afterthreshold(Pos_th)=0;
            
            Q=8;% 1 byte (value representation)
            Max_SigAfter=max(Sig_Afterthreshold);
            Min_SigAfter=min(Sig_Afterthreshold);
            Sig_Quant=round((-1+2^Q)*(Sig_Afterthreshold-Min_SigAfter)/(Max_SigAfter-Min_SigAfter));
            
            Sig_RLE=rle(Sig_Quant)
            
            Sig_Irle=irle(Sig_RLE);% Inverse Run Length Codising 
            Sig_IQuant=((Max_SigAfter-Min_SigAfter)/(-1+2^Q))*Sig_Irle+Min_SigAfter;% Inverse Quantization
            Rec_Sig=ifwht(Sig_IQuant);%Décompressed Signal
            save('Rec_Sig');
            caxes(handles.axes2); 
            
            plot(Rec_Sig);
        case 1 %Lossless
            wht_Sig=fwht(signal);
            
            Sig_RLE=rle(wht_Sig);
            Sig_RLE=single(Sig_RLE);
            Sig_RLE=double(Sig_RLE);
            
            Sig_Irle=irle(Sig_RLE);% Inverse Run Length Codising 
            Rec_Sig=ifwht(Sig_Irle);%Décompressed Signal
            save('Rec_Sig');
            axes(handles.axes2); 
            plot(Rec_Sig);
    end
case 3 %W2
    switch Lflag
        case 0 %Lossy
             y=1:1:length(signal);

            DEC = mdwtdec('c',signal,5,'sym4');

            [XC,DECCMP,THRESH] = mswcmp('cmp',DEC,'bal_sn');
            mu = 0;
            sigma = 1;
            pd = makedist('Normal','mu',mu,'sigma',sigma);
            p = pdf (pd,XC);
                p1=p/sum(p);
            XCu = unique(XC);
            len=1:1:length(XCu);
            prob = 1:length(len);
            counter = 0;
            for i = 1:length(len)
                for j = 1:length(y)

                        if (XCu(i)==XC(j))
                          counter = counter+1;
                        end         
                end
                prob(i)=counter/length(y);
                counter = 0;
            end
            [dict,~] = huffmandict(XCu,prob);
            enc = huffmanenco(XC,dict);
            dsig = huffmandeco(enc,dict)
            save('dsig')
  
            axes(handles.axes2);
            plot ([dsig;signal-dsig]');  


            legend('signal ','Diff');
        case 1 %Lossless
    end
end


% --- Executes on selection change in popMenu.
function popMenu_Callback(hObject, eventdata, handles)

contents = cellstr(get(hObject,'String'))% returns popMenu contents as cell array
global Lflag
Lselection = contents{get(hObject,'Value')} %returns selected item from popMenu
if (strcmp(Lselection,'Lossy'))
Lflag = 0;
elseif (strcmp(Lselection,'Lossless'))
Lflag = 1;
end
    % --- Executes during object creation, after setting all properties.
function popMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popMenuBF.
function popMenuBF_Callback(hObject, eventdata, handles)
% hObject    handle to popMenuBF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String')) %returns popMenuBF contents as cell array
global BFflag
BFselection = contents{get(hObject,'Value')} %returns selected item from popMenuBF
if (strcmp(BFselection,'DCT'))
BFflag = 0;
elseif (strcmp(BFselection,'DFT'))
BFflag = 1;
elseif (strcmp(BFselection,'W1'))
BFflag = 2;
elseif (strcmp(BFselection,'W2'))
BFflag = 3;
end

% --- Executes during object creation, after setting all properties.
function popMenuBF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popMenuBF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
global signal
clc
  [filename pathname] = uigetfile({'*.mat;*.xlsx;*.wav';'*.*'},'Select File');
   handles.fullpathname = strcat(pathname, filename);
   [filepath,name,ext]=fileparts(filename)
   if(strcmp(ext,'.wav'))
       audioread(filename);
       signal=ans;
   elseif (strcmp(ext,'.mat'))
     load(filename);
        if (strcmp(filename,'ult_sig.mat'))
        signal = ult_sig ;
        else
            signal = val;
        end
   else
      xlsread(filename); 
      signal= ans ;
   end 
   axes(handles.axes1);
   plot(signal);