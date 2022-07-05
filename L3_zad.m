%% WCZYTANIE FILMU

N = 300; % liczba ramek do wczytania (pierwsze 10 sekund nagrania)
br = zeros(1, N); % wektor jasności

FPS = 30;   % częstotliwość próbkowania [Hz]
T = 1/FPS;  % okres próbkowania [s]
t = (0:N-1)*T; % oś czasu

v = VideoReader(['pulsKamila.mp4']); % załadowanie i odczytanie pliku video

%% UZYSKANIE WYKRESU JASNOŚCI

% wczytanie pierwszych N obrazów i analiza jasności
for i=1:N
    I = rgb2gray(read(v,i)); % wczytanie klatki w skali szarości
    br(i) = mean(I, 'all'); % średnia jasność całego obrazu w danej klatce
end

br = br - mean(br); % przesunięcie o stałą

figure % wykres czasowy bez filtracji
plot(t, br)
xlabel("czas [s]")
ylabel("jasność nagrania")
title("Jasność nagrania w czasie (bez filtrów)")
print("jasność_bez.png","-dpng","-r280")

%% FILTROWANIE WYNIKU

N_filter = 5;

% Filtr uśredniający:
f_mean = ones(1,N_filter)/N_filter;
c_mean = conv(br, f_mean, "same");

figure % wykres czasowy
plot(t, c_mean)
xlabel("czas [s]")
ylabel("jasność nagrania")
title(strcat("Jasność nagrania w czasie (filtr uśredniający, N=", string(N_filter), ")"))
print("jasność_z_usrednieniem.png","-dpng","-r280")

% Filtr Gaussa:
f_gauss = fspecial('gaussian', [1, N_filter], 3);
c_gauss = conv(br, f_gauss, "same");

figure % wykres czasowy
plot(t, c_gauss)
xlabel("czas [s]")
ylabel("jasność nagrania")
title(strcat("Jasność nagrania w czasie (filtr Gaussa, N=", string(N_filter), ")"))
print("jasność_z_gaussem.png","-dpng","-r280")

% Detekcja zbocza:
f_cliff = linspace(1, -1, N_filter);
c_cliff = conv(br, f_cliff, "same");

figure % wykres czasowy
plot(t, c_cliff)
xlabel("czas [s]")
ylabel("jasność nagrania")
title("Jasność nagrania w czasie (detekcja zbocza)")
print("jasność_zbocza.png","-dpng","-r280")

% Autokorelacja
[r, lags] = xcorr(br);
r = r(lags >= 0);
lags = lags(lags>=0);

figure % wykres czasowy
plot(lags, r)
xlabel("czas [s]")
ylabel("jasność nagrania")
title("Jasność nagrania w czasie (autokorelacja)")
print("jasność_autokorelacja.png","-dpng","-r280")


%% TRANSFORMATA FOURIERA 

Y = fft(br);    % szybka transformata Fouriera
A = abs(Y);     % amplituda sygnału
A = A/N;        % normalizacja amplitudy
A = A(1:N/2+1); % wycięcie istotnej części spektrum
A(2:end-1) = 2*A(2:end-1);

Y_mean = fft(c_mean);    % szybka transformata Fouriera
A_mean = abs(Y_mean);     % amplituda sygnału
A_mean = A_mean/N;        % normalizacja amplitudy
A_mean = A_mean(1:N/2+1); % wycięcie istotnej części spektrum
A_mean(2:end-1) = 2*A_mean(2:end-1);

Y_gauss = fft(c_gauss);    % szybka transformata Fouriera
A_gauss = abs(Y_gauss);     % amplituda sygnału
A_gauss = A_gauss/N;        % normalizacja amplitudy
A_gauss = A_gauss(1:N/2+1); % wycięcie istotnej części spektrum
A_gauss(2:end-1) = 2*A_gauss(2:end-1);

Y_cliff = fft(c_cliff);    % szybka transformata Fouriera
A_cliff = abs(Y_cliff);     % amplituda sygnału
A_cliff = A_cliff/N;        % normalizacja amplitudy
A_cliff = A_cliff(1:N/2+1); % wycięcie istotnej części spektrum
A_cliff(2:end-1) = 2*A_cliff(2:end-1);

Y_r = fft(r);    % szybka transformata Fouriera
A_r = abs(Y_r);     % amplituda sygnału
A_r = A_r/N;        % normalizacja amplitudy
A_r = A_r(1:N/2+1); % wycięcie istotnej części spektrum
A_r(2:end-1) = 2*A_r(2:end-1);

f_step = FPS/N;     % zmiana częstotliwości
f = 0:f_step:FPS/2; % oś częstotliwości

figure % wykres amplitudowy
plot(f, A)
xlabel("częstotliwość [Hz]")
ylabel("amplituda")
title("Wyniki FFT - bez filtrów")
print("FFT_bez.png","-dpng","-r280")

figure % wykres amplitudowy
plot(f, A_mean)
xlabel("częstotliwość [Hz]")
ylabel("amplituda")
title("Wyniki FFT - filtr uśredniający")
print("FFT_mean.png","-dpng","-r280")

figure % wykres amplitudowy
plot(f, A_gauss)
xlabel("częstotliwość [Hz]")
ylabel("amplituda")
title("Wyniki FFT - filtr Gaussa")
print("FFT_gauss.png","-dpng","-r280")

figure % wykres amplitudowy
plot(f, A_cliff)
xlabel("częstotliwość [Hz]")
ylabel("amplituda")
title("Wyniki FFT - detekcja zbocza")
print("FFT_zbocze.png","-dpng","-r280")

figure % wykres amplitudowy
plot(f, A_r)
xlabel("częstotliwość [Hz]")
ylabel("amplituda")
title("Wyniki FFT - autokorelacja")
print("FFT_autokorelacja.png","-dpng","-r280")

%% ODCZYTANIE CZĘSTOTLIWOŚCI

[max_shade_change, heart_rate_sample] = maxk(A, 1); % znalezienie indeksu dla tętna serca
heart_BPM = 60*f(heart_rate_sample) % przeliczenie na BPM

[max_shade_change, heart_rate_sample] = maxk(A_mean, 1); % znalezienie indeksu dla tętna serca
heart_BPM_mean = 60*f(heart_rate_sample) % przeliczenie na BPM

[max_shade_change, heart_rate_sample] = maxk(A_gauss, 1); % znalezienie indeksu dla tętna serca
heart_BPM_gauss = 60*f(heart_rate_sample) % przeliczenie na BPM

[max_shade_change, heart_rate_sample] = maxk(A_cliff, 1); % znalezienie indeksu dla tętna serca
heart_BPM_cliff = 60*f(heart_rate_sample) % przeliczenie na BPM

[max_shade_change, heart_rate_sample] = maxk(A_r, 1); % znalezienie indeksu dla tętna serca
heart_BPM_autor = 60*f(heart_rate_sample) % przeliczenie na BPM

%% WYKRESY ŁĄCZONE

figure
subplot(2,3,1)
plot(t, br)
title("bez filtrów")
subplot(2,3,2)
plot(t, c_mean)
title("filtr uśredniający")
subplot(2,3,3)
plot(t, c_gauss)
title("filtr Gaussa")
subplot(2,3,5)
plot(t, c_cliff)
title("detekcja zbocza")
subplot(2,3,6)
plot(t, r)
title("autokorelacja")
print("czasówki_razem.png","-dpng","-r280")

figure
subplot(2,3,1)
plot(f, A)
title("FFT bez filtrów")
subplot(2,3,2)
plot(f, A_mean)
title("FFT filtr uśredniający")
subplot(2,3,3)
plot(f, A_gauss)
title("FFT filtr Gaussa")
subplot(2,3,5)
plot(f, A_cliff)
title("FFT detekcja zbocza")
subplot(2,3,6)
plot(f, A_r)
title("FFT autokorelacja")
print("FFT_razem.png","-dpng","-r280")


