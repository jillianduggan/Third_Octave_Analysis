# inputs:
#   1. wav file (in apostrophes)
#   2. hydrophone sensitivity in dB
#   3. start time in seconds
#   4. stop time (if you want it to be the end of the wav file, write [])

# Other things to note:
#   - The fft window size is set to 65536 samples
#   - Window overlap is 50%
#   - If file is shorter than 65536 samples, then window is set to the total number of samples
#   - The samples that do not fit into the last window are left out of the analysis (this is what octave does too, in calculating the power spectrum)



function third_oct(file,sens,start,stop)

disp(' ');

# Send error if first input is not a string
if typeinfo(file) != 'sq_string'
error('Wav file must be a string. Enclose text in apostrophes or quotes.')
end


# Read wav file
[y,fs] = audioread(file);
n = length(y);
t = num2str(n/fs);

# Send error if wav file is shorter than start time
if start >= n/fs
error(["Recording is shorter than the chosen start time.\n       Recording length is " t " seconds."])
end

# Calculate stop time if stop is []
if stop == 'end'
stop = n/fs;
end


# Send error if start time >= stop time
 if start >= stop
error('Stop time must be greater than the start time.')
end


# Take the portion of wav file between the start and stop times
y = y(start*fs+1:stop*fs);



# new length
n = length(y);

# FFT window size
w = 65536;

# Number of windows
nwin = floor(n/(w/2))-1;
if nwin == 0
w = length(y);
nwin = 1;
s = num2str(w);
disp(['Warning: wav file shorter than fft size 65536. New fft size = ' s '.'])
end

# Divide signal into windows and apply Hanning window to each
Y = zeros(nwin,w);
for i = 1:nwin
Y(i,:) = hanning(w).*y((i-1)*w/2+1:(i+1)*w/2);
end

# Fourier transform each window
C = zeros(nwin,w);
for i = 1:nwin
C(i,:) = abs(fft(Y(i,:),w,2)).^2;
end

# Average over windows
S = zeros(1,w);
S = sum(C,1)/nwin;

# Normalize the FFT
S = S/(sqrt(2)*w);

# Keep left half of the spectrum
S = S(1:floor(w/2));

# New length of spectrum
N = w/2;

# Power spectrum in dB
Ps = 20*log10(S);

# Adjust for sensitivity
S = S*10^(-sens/20);

# Density
S = S/sqrt(fs/w);


# Spectrum frequency bins
F = (1:w/2)*fs/w;


# Standard third octave values
# (goes up to 250 kHz)
# the i^th band is (T(i-1), T(i))
T = [0;
     16;
     20;
     25;
     32;
     40;
     50;
     63;
     80;
     100;
     125;
     160;
     200;
     250;
     315;
     400;
     500;
     630;
     800;
     1000;
     1250;
     1600;
     2000;
     2500;
     3150;
     4000;
     5000;
     6300;
     8000;
     10000;
     12500;
     16000;
     20000;
     25000;
     31500;
     40000;
     50000;
     63000;
     80000;
     100000;
     125000;
     160000;
     200000;
     250000];



# Calculated third octave bands, base 10
# This is the standard way...
# ...the values in T correspond to rounded centre frequencies
# fcentre = actual band centre frequeny
# fupper = band upper limit
# flower = band lower limit
fcentre = 10.^(0.1.*[12:55]);
fd = 10^0.05;
fupper = fcentre * fd;
flower = fcentre / fd;
flower(1) = 1;


    
# Remove unused bins
J = fcentre < F(N);
fcentre = fcentre(J);
m = length(fcentre);
T = T(1:m);
      

      
# sum spectrum values over 1/3 octave bins
      Z = zeros(m,1);
      for i = 1:m
      J = (flower(i)<= F) & (F < fupper(i));
      Z(i) = sum(S(J));
      end
      
    

# Convert to dB
      for i = 1:m
      if Z(i) != 0
      Z(i) = 20*log10(Z(i));
      else
      Z(i) = 0;
      end
      end


# Prepare for writing to text file
      T = T(1:m);
      TOA = [T Z];
      a = index(file,'.wav');
      
# To suppress file info, add '-ascii' as the first argument
      save(['third_octave_' sprintf(strtrunc(file,a-1))],'TOA')

disp(' ');
disp(['A 1/3 octave analysis of ' sprintf(file) ' was performed and saved in the text file third_octave_' sprintf(strtrunc(file,a-1)) '.'])
disp(' ');
      
# figure
# scatter(1:length(Z),Z)


