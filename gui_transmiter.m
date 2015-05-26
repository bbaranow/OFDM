function varargout = gui_transmiter(varargin)
% GUI_TRANSMITER MATLAB code for gui_transmiter.fig
%      GUI_TRANSMITER, by itself, creates a new GUI_TRANSMITER or raises the existing
%      singleton*.
%
%      H = GUI_TRANSMITER returns the handle to a new GUI_TRANSMITER or the handle to
%      the existing singleton*.
%
%      GUI_TRANSMITER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TRANSMITER.M with the given input arguments.
%
%      GUI_TRANSMITER('Property','Value',...) creates a new GUI_TRANSMITER or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_transmiter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_transmiter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_transmiter

% Last Modified by GUIDE v2.5 26-May-2015 20:30:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_transmiter_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_transmiter_OutputFcn, ...
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

% --- Executes just before gui_transmiter is made visible.
function gui_transmiter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_transmiter (see VARARGIN)

% Choose default command line output for gui_transmiter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes gui_transmiter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_transmiter_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function carrier_count_CreateFcn(hObject, eventdata, handles)
% hObject    handle to carrier_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function carrier_count_Callback(hObject, eventdata, handles)
% hObject    handle to carrier_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of carrier_count as text
%        str2double(get(hObject,'String')) returns contents of carrier_count as a double
carrier_count = str2double(get(hObject, 'String'));
if isnan(carrier_count)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new carrier_count value
handles.metricdata.carrier_count = carrier_count;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function bits_per_symbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bits_per_symbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bits_per_symbol_Callback(hObject, eventdata, handles)
% hObject    handle to bits_per_symbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bits_per_symbol as text
%        str2double(get(hObject,'String')) returns contents of bits_per_symbol as a double
bits_per_symbol = str2double(get(hObject, 'String'));
if isnan(bits_per_symbol)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new bits_per_symbol value
handles.metricdata.bits_per_symbol = bits_per_symbol;
guidata(hObject,handles)

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IFFT_bin_length=256;
carrier_count=handles.metricdata.carrier_count;
symbols_per_carrier=handles.metricdata.symbols_per_carrier;
bits_per_symbol=handles.metricdata.bits_per_symbol;
baseband_out_length =  carrier_count * symbols_per_carrier * bits_per_symbol; 
carriers = (1:carrier_count) + (floor(IFFT_bin_length/4) - floor(carrier_count/2));
conjugate_carriers = IFFT_bin_length - carriers + 2;
%
%
%--------1---------2---------3---------4---------5---------6---------7---------8
%
% TRANSMIT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
%
% Generate a random binary output signal:
%   - a row of uniform random numbers (between 0 and 1), rounded to 0 or 1
%   - this will be the baseband signal which is to be transmitted.
%
%baseband_out = round(rand(1,baseband_out_length));
baseband_out = importdata('danewe.txt');
%
% Convert to 'modulo N' integers where N = 2^bits_per_symbol
%   - this defines how many states each symbol can represent
%   - first, make a matrix with each column representing consecutive bits 
%     from the input stream and the number of bits in a column equal to the
%     number of bits per symbol
%   - then, for each column, multiply each row value by the power of 2 that 
%     it represents and add all the rows
%   - for example:  input 0 1 1 0 0 0 1 1 1 0
%                   bits_per_symbol = 2
%                   convert_matrix = 0 1 0 1 1
%                                    1 0 0 1 0
%
%                   modulo_baseband = 1 2 0 3 2
%
convert_matrix = reshape(baseband_out, bits_per_symbol, length(baseband_out)/bits_per_symbol);
for k = 1:(length(baseband_out)/bits_per_symbol)
    modulo_baseband(k) = 0;
    for i = 1:bits_per_symbol
        modulo_baseband(k) = modulo_baseband(k) + convert_matrix(i,k)*2^(bits_per_symbol-i);
    end
end
%
%--------1---------2---------3---------4---------5---------6---------7---------8
%
% Serial to Parallel Conversion
%   - convert the serial modulo N stream into a matrix where each column 
%     represents a carrier and each row represents a symbol
%   - for example:
%
%      serial input stream = a b c d e f g h i j k l m n o p
%
%      parallel carrier distribution =
%                           C1/s1=a  C2/s1=b  C3/s1=c C4/s1=d
%                           C1/s2=e  C2/s2=f  C3/s2=g C4/s2=h
%                           C1/s3=i  C2/s3=j  C3/s3=k C4/s3=l
%                              .        .        .       .
%                              .        .        .       .
%
carrier_matrix = reshape(modulo_baseband, carrier_count, symbols_per_carrier)';
%
% Apply differential coding to each carrier string
%   - append an arbitrary start symbol (let it be 0, that works for all 
%     values of bits_per_symbol) (note that this is done using a vertical
%     concatenation [x;y] of a row of zeros with the carrier matrix, sweet!)
%   - perform modulo N addition between symbol(n) and symbol(n-1) to get the 
%     coded value of symbol(n)
%   - for example:
%                   bits_per_symbol = 2 (modulo 4)
%                   symbol stream =  3 2 1 0 2 3
%                   start symbol = 0
%
%                   coded symbols = 0 + 3 = 3
%                                   3 + 2 = 11 = 1
%                                   1 + 1 = 2
%                                   2 + 0 = 2
%                                   2 + 2 = 10 = 0
%                                   0 + 3 = 3
%
%                   coded stream =  0 3 1 2 2 0 3
%
%--------1---------2---------3---------4---------5---------6---------7---------8
%
carrier_matrix = [zeros(1,carrier_count);carrier_matrix];
for i = 2:(symbols_per_carrier + 1)
    carrier_matrix(i,:) = rem(carrier_matrix(i,:)+carrier_matrix(i-1,:),2^bits_per_symbol);
end
%
% Convert the differential coding into a phase
%   - each phase represents a different state of the symbol
%   - for example:
%                 bits_per_symbol = 2 (modulo 4)
%                 symbols = 0 3 2 1
%                 phases =
%                           0 * 2pi/4 = 0 (0 degrees)
%                           3 * 2pi/4 = 3pi/2 (270 degrees)
%                           2 * 2pi/4 = pi (180 degrees)
%                           1 * 2pi/4 = pi/2 (90 degrees)
%
carrier_matrix = carrier_matrix * ((2*pi)/(2^bits_per_symbol));
%
% Convert the phase to a complex number
%   - each symbol is given a magnitude of 1 to go along with its phase 
%     (via the ones(r,c) function)
%   - it is then converted from polar to cartesian (complex) form
%   - the result is 2 matrices, X with the real values and Y with the imaginary
%   - each X column has all the real values for a carrier, and each Y column 
%     has the imaginary values for a carrier
%   - a single complex matrix is then generated taking X for real and 
%     Y for imaginary
%
[X,Y] = pol2cart(carrier_matrix, ones(size(carrier_matrix,1),size(carrier_matrix,2)));
complex_carrier_matrix = complex(X,Y);
%
%--------1---------2---------3---------4---------5---------6---------7---------8
%
% Assign each carrier to its IFFT bin
%   - each row of complex_carrier_matrix represents one symbol period, with
%     a symbol for each carrier
%   - a matrix is generated to represent all IFFT bins (columns) and all 
%     symbols (rows)
%   - the phase modulation for each carrier is then assigned to the 
%     appropriate bin
%   - the conjugate of the phase modulation is then assigned to the 
%     appropriate bin
%      - the phase modulation bins and their conjugates are symmetric about 
%        the Nyquist frequency in the IFFT bins
%      - since the first bin is DC, the Nyquist Frequency is located 
%        at (number of bins/2) + 1
%      - symmetric conjugates are generated so that when the signal is 
%        transformed to the time domain, the time signal will be real-valued
%   - example
%      - 1024 IFFT bins
%      - bin 513 is the center (symmetry point)
%      - bin 1 is DC
%      - bin 514 is the complex conjugate of bin 512
%      - bin 515 is the complex conjugate of bin 511
%      - ....
%      - bin 1024 is the complex conjugate of bin 2 (if all bins 
%        were used as carriers)
%      - So, bins 2-512 map to bins 1024-514
%
IFFT_modulation = zeros(symbols_per_carrier + 1, IFFT_bin_length);
IFFT_modulation(:,carriers) = complex_carrier_matrix;
IFFT_modulation(:,conjugate_carriers) = conj(complex_carrier_matrix);
%
% PLOT BASIC FREQUENCY DOMAIN REPRESENTATION
%
%--------1---------2---------3---------4---------5---------6---------7---------8
%
%figure (1)
%stem(0:IFFT_bin_length-1, abs(IFFT_modulation(2,1:IFFT_bin_length)),'b*-')
%grid on
%axis ([0 IFFT_bin_length -0.5 1.5])
%ylabel('Magnitude')
%xlabel('IFFT Bin')
%title('OFDM Carrier Frequency Magnitude')
%figure (2)
%plot(0:IFFT_bin_length-1, (180/pi)*angle(IFFT_modulation(2,1:IFFT_bin_length)), 'go')
%hold on
%stem(carriers-1, (180/pi)*angle(IFFT_modulation(2,carriers)),'b*-')
%stem(conjugate_carriers-1, (180/pi)*angle(IFFT_modulation(2,conjugate_carriers)),'b*-')
%axis ([0 IFFT_bin_length -200 +200])
%grid on
%ylabel('Phase (degrees)')
%xlabel('IFFT Bin')
%title('OFDM Carrier Phase')
% END OF PLOTTING
%
% Transform each period's spectrum (represented by a row of carriers) to the 
% time domain via IFFT
%
time_wave_matrix = ifft(IFFT_modulation');
time_wave_matrix = time_wave_matrix';
%
% PLOT OFDM SIGNAL FOR ONE SYMBOL PERIOD
%   - first plot is direct IFFT of first OFDM symbol (index 2)
%      - index 1 is the 'start' symbol due to differential implementation
%   - second plot strips out each carrier and plots them on the
%     same graph
%      - only works for carrier count <= 16 due to colors variable (more
%        than 16 would really be legible anyway
%
%figure (3)
%plot(0:IFFT_bin_length-1,time_wave_matrix(2,:))
%grid on
%ylabel('Amplitude')
%xlabel('Time')
%title('OFDM Time Signal, One Symbol Period')
%
%colors = ['r' 'g' 'b' 'k' 'r' 'g' 'b' 'k' 'r' 'g' 'b' 'k' 'r' 'g' 'b' 'k'];
%for f = 1:carrier_count
%    temp_bins(1:IFFT_bin_length)=0+0j;
%    temp_bins(carriers(f))=IFFT_modulation(2,carriers(f));
%    temp_bins(conjugate_carriers(f))=IFFT_modulation(2,conjugate_carriers(f));
%    temp_time = ifft(temp_bins');
%    figure(4)
%    plot(0:IFFT_bin_length-1, temp_time, colors(f))
%   hold on
%end
%grid on
%ylabel('Amplitude')
%xlabel('Time')
%title('Separated Time Waveforms Carriers')
%
% END OF PLOTTING
%
%--------1---------2---------3---------4---------5---------6---------7---------8
%
% Apply a Window Function to each time waveform
%   - NOTE THAT WINDOWING IS CURRENTLY COMMENTED OUT, i.e. NO WINDOWING
%   - each time waveform (row of time_wave_matrix) represents one symbol 
%     period for all carriers
%   - the IFFT result has discontinuities at each end 
%   - when the time waveforms are serialized (concatenated), the discontinuites
%     will introduce unwanted frequency components
%   - the window function deemphasizes the signal at the end 
%     points (at the discontinuites)
%      - this reduces the effects of the discontinuities
%      - it also distorts the desired frequency response (undesired side effect)
%   - between Blackman, Hanning, and Hamming:  Hamming introduces less distortion
%   - note that the transpose of the Hamming function is 
%     used (because a row vector is needed)
%
% Since all imaginary values of time_wave_matrix are practically equal to zero,
% only the real part is retained for windowing.
%
for i = 1:symbols_per_carrier + 1
    %windowed_time_wave_matrix(i,:) = real(time_wave_matrix(i,:)) .* hamming(IFFT_bin_length)';
    windowed_time_wave_matrix(i,:) = real(time_wave_matrix(i,:));
end
%
% Serialize the modulating waveform
%   - sequentially take each row of windowed_time_wave_matrix and construct a row vector
%   - the row vector will be the modulating signal
%   - note that windowed_time_wave_matrix is transposed, this is to account for the way the
%     Matlab 'reshape' function works (reshape takes the columns of the target matrix and 
%     appends them sequentially)
% 
ofdm_modulation = reshape(windowed_time_wave_matrix', 1, IFFT_bin_length*(symbols_per_carrier+1));
%
% PLOT OFDM SIGNAL (time)
%
%temp_time = IFFT_bin_length*(symbols_per_carrier+1);
%figure (5)
%plot(0:temp_time-1,ofdm_modulation)
%grid on
%ylabel('Amplitude (volts)')
%xlabel('Time (samples)')
%title('OFDM Time Signal')
%
% PLOT OFDM SIGNAL (spectrum)
%symbols_per_average = ceil(symbols_per_carrier/5);
%avg_temp_time = IFFT_bin_length*symbols_per_average;
%averages = floor(temp_time/avg_temp_time);
%average_fft(1:avg_temp_time) = 0;
%for a = 0:(averages-1)
%    subset_ofdm = ofdm_modulation(((a*avg_temp_time)+1):((a+1)*avg_temp_time));
%    subset_ofdm_f = abs(fft(subset_ofdm));
%    average_fft = average_fft + (subset_ofdm_f/averages);
%end
%average_fft_log = 20*log10(average_fft);
%figure (6)
%plot((0:(avg_temp_time-1))/avg_temp_time, average_fft_log)
%hold on
%plot(0:1/IFFT_bin_length:1, -35, 'rd')
%grid on
%axis([0 0.5 -40 max(average_fft_log)])
%ylabel('Magnitude (dB)')
%xlabel('Normalized Frequency (0.5 = fs/2)')
%title('OFDM Signal Spectrum')
%
% ENDPLOT
%
%--------1---------2---------3---------4---------5---------6---------7---------8
%
% Upconversion to RF 
%
% For this model, the baseband will be inserted directly into the channel
% without conversion to RF frequencies.
%
Tx_data = ofdm_modulation;
dlmwrite('wynik.txt', Tx_data);
% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IFFT_bin_length=256;
carrier_count=handles.metricdata.carrier_count;
symbols_per_carrier=handles.metricdata.symbols_per_carrier;
bits_per_symbol=handles.metricdata.bits_per_symbol;
baseband_out_length =  carrier_count * symbols_per_carrier * bits_per_symbol; 
carriers = (1:carrier_count) + (floor(IFFT_bin_length/4) - floor(carrier_count/2));
conjugate_carriers = IFFT_bin_length - carriers + 2;

Tx_data = importdata('wynik.txt');
Tx_signal_power = var(Tx_data);
%

Rx_Data = Tx_data;
%
%
% RECEIVE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%
%
% Convert the serial input data stream to parallel (according to symbol length 
% and number of symbols) 
%   - each column is a symbol period
%   - the length of each symbol (samples per symbol) is the length of the 
%     IFFT that was used to generate it
%
Rx_Data_matrix = reshape(Rx_Data, IFFT_bin_length, symbols_per_carrier + 1);
%
% Transform each symbol from time to frequency domain
%   - take the fft of each column
%
Rx_spectrum = fft(Rx_Data_matrix);%
% PLOT BASIC FREQUENCY DOMAIN REPRESENTATION
%
%--------1---------2---------3---------4---------5---------6---------7---------8
%
%figure (7)
%stem(0:IFFT_bin_length-1, abs(Rx_spectrum(1:IFFT_bin_length,2)),'b*-')
%grid on
%axis ([0 IFFT_bin_length -0.5 1.5])
%ylabel('Magnitude')
%xlabel('FFT Bin')
%title('OFDM Receive Spectrum, Magnitude')
%figure (8)
%plot(0:IFFT_bin_length-1, (180/pi)*angle(Rx_spectrum(1:IFFT_bin_length,2)), 'go')
%hold on
%stem(carriers-1, (180/pi)*angle(Rx_spectrum(carriers,2)),'b*-')
%stem(conjugate_carriers-1, (180/pi)*angle(Rx_spectrum(conjugate_carriers,2)),'b*-')
%axis ([0 IFFT_bin_length -200 +200])
%grid on
%ylabel('Phase (degrees)')
%xlabel('FFT Bin')
%title('OFDM Receive Spectrum,  Phase')
%
% END OF PLOTTING


%--------1---------2---------3---------4---------5---------6---------7---------8
%
% Extract the carrier FFT bins
%   - only keep the fft bins that are used as carriers
%   - take the transpose of the result so that each column will represent 
%     a carrier
%      - this is in preparation for using the diff( ) function later to decode 
%        differential encoding
%      - format following this operation is:
%
%               C1-s1  C2-s1  C3-s1 ...
%               C1-s2  C2-s2  C3-s2 ...
%               C1-s3  C2-s3  C3-s3 ...
%                   .      .      .
%                   .      .      .
%         
%   - IMPORTANT MATLAB NOTE CONCERNING TRANSPOSING AND CONJUGATION
%     - it appears that each time a matrix is transposed, the conjugate of 
%       each value is taken
%     - if an even number of transposes are done, then it is transparent
%     - obviously, this does not affect real numbers
%
Rx_carriers = Rx_spectrum(carriers,:)';
%
%--------1---------2---------3---------4---------5---------6---------7---------8
%
% PLOT EACH RECEIVED SYMBOL
%
%figure (9)
Rx_phase_P = angle(Rx_carriers);
Rx_mag_P = abs(Rx_carriers);
%polar(Rx_phase_P, Rx_mag_P,'bd');
%
% END PLOT
%
% Find the phase (angle) of each FFT bin (each carrier)
%   - convert from radians to degrees
%   - normalize phase to be between 0 and 359 degrees
%
Rx_phase = angle(Rx_carriers)*(180/pi);
phase_negative = find(Rx_phase < 0);
Rx_phase(phase_negative) = rem(Rx_phase(phase_negative)+360,360);
%
% Extract phase differences (from the differential encoding)
%   - the matlab diff( ) function is perfect for this operation
%   - again, normalize the result to be between 0 and 359 degrees
%
Rx_decoded_phase = diff(Rx_phase);
phase_negative = find(Rx_decoded_phase < 0);
Rx_decoded_phase(phase_negative) = rem(Rx_decoded_phase(phase_negative)+360,360);
%
%--------1---------2---------3---------4---------5---------6---------7---------8
%
% Convert phase to symbol
%   - calculate the base phase which is the phase difference between each 
%     consecutive symbol
%   - for example, if there are 2 bits per symbol, base phase is 90 and the 
%     symbols are represented by 0, 90, 180, and 270 degrees
%   - calculate the maximum deviation from the base phase that will still be 
%     decoded as base phase
%   - for example, if base phase is 90, then delta phase is 45, and anything 
%     within 45 degrees of base phase is accepted as base phase
%      - continuing the above example, a symbol represented by 180 will be 
%        decoded as 180 as long as it is within the range 180+45 and 180-45
%   - generate a symbol matrix where the results of the phase decode 
%     will be placed
%      - note that since the matrix is created as a zero matrix, then the zero
%        values do not have to be decoded
%      - zero is therefore the default value after decoding, if a value is not 
%        decoded as anything else, then it is zero
%      - this actually save alot of trouble since the zero phase spans the 
%        lowest and highest phase values and therefore requires special 
%        processing
%      - it is also efficient in that it eliminates a pass through the loop
%
base_phase = 360/2^bits_per_symbol;
delta_phase = base_phase/2;
Rx_decoded_symbols = zeros(size(Rx_decoded_phase,1),size(Rx_decoded_phase,2));
%
for i = 1:(2^bits_per_symbol - 1)  
    center_phase = base_phase*i;
    plus_delta = center_phase+delta_phase;
    minus_delta = center_phase-delta_phase;
    decoded = find((Rx_decoded_phase <= plus_delta) & (Rx_decoded_phase > minus_delta));
    Rx_decoded_symbols(decoded)=i;
end
%
% Convert the matrix into a serial symbol stream
%
Rx_serial_symbols = reshape(Rx_decoded_symbols',1,size(Rx_decoded_symbols,1)*size(Rx_decoded_symbols,2));
%
% Convert the symbols to binary
%
for i = bits_per_symbol: -1: 1
    if i ~= 1
        Rx_binary_matrix(i,:) = rem(Rx_serial_symbols,2);
        Rx_serial_symbols = floor(Rx_serial_symbols/2);
    else
        Rx_binary_matrix(i,:) = Rx_serial_symbols;
    end
end
baseband_in = reshape(Rx_binary_matrix,1,size(Rx_binary_matrix,1)*size(Rx_binary_matrix,2));
dlmwrite('danewy.txt',baseband_in);
%
% Find bit errors
%
%--------1---------2---------3---------4---------5---------6---------7---------8

% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (hObject == handles.english)
    set(handles.text4, 'String', 'lb/cu.in');
    set(handles.text5, 'String', 'cu.in');
    set(handles.text6, 'String', 'lb');
else
    set(handles.text4, 'String', 'kg/cu.m');
    set(handles.text5, 'String', 'cu.m');
    set(handles.text6, 'String', 'kg');
end

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

handles.metricdata.carrier_count = 0;
handles.metricdata.bits_per_symbol  = 0;
handles.metricdata.symbols_per_carrier  = 0;
set(handles.carrier_count, 'String', handles.metricdata.carrier_count);
set(handles.bits_per_symbol,  'String', handles.metricdata.bits_per_symbol);
%set(handles.symbols_per_carrier, 'String',handles.metricdata.symbols_per_carrier);

% Update handles structure
guidata(handles.figure1, handles);



function symbols_per_carrier_Callback(hObject, eventdata, handles)
% hObject    handle to symbols_per_carrier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of symbols_per_carrier as text
%        str2double(get(hObject,'String')) returns contents of symbols_per_carrier as a double
symbols_per_carrier = str2double(get(hObject, 'String'));
if isnan(symbols_per_carrier)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new bits_per_symbol value
handles.metricdata.symbols_per_carrier = symbols_per_carrier;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function symbols_per_carrier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to symbols_per_carrier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over calculate.
function calculate_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
