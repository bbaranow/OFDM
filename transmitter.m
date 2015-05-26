%--------1---------2---------3---------4---------5---------6---------7---------8
% OFDM Simulation
%
%    - ECE 5664 Project
%    - Erich Cosby
%    - December 4, 2001
%    - Virginia Tech, Nothern Va. Center
%
% Based on the matlab work of Eric Lawrey in his 1997 BSEE thesis,
% "The suitability of OFDM as a modulation technique for wireless 
%  telecommunications, with a CDMA comparison", 
% refer to www.eng.jcu.edu.au/eric/thesis/Thesis.htm.
%
%
clear all;
close all;
%

%
%if defaults == 1
    IFFT_bin_length = 1024;     % IFFT bin count for Tx, FFT bin count for Rx
    carrier_count = 200;        % number of carriers
    bits_per_symbol = 4;        % bits per symbol
    symbols_per_carrier = 50;   % symbols per carrier
    SNR = 10;                   % channel signal to noise ratio (dB)                           
%else
%    IFFT_bin_length = input('IFFT bin length = ');
%    carrier_count = input('carrier count = ');
%    bits_per_symbol = input('bits per symbol = ');
%    symbols_per_carrier = input('symbols per carrier =');
%    SNR = input('SNR = ');
%end
%
% Derived parameters
%
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
baseband_out = round(rand(1,baseband_out_length));
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