# Projekt poprawkowy na przedmiot Elektronika
# Autor: Jakub Styn kod: 330
# magiczna funkcja IPythona umożliwiająca kreślenie w wierszu, w którym wykresy/wykresy będą wyświetlane tuż pod komórką,
# w której zapisano polecenia kreślenia.
%matplotlib inline
import math  # importuję bliotekę math zawierającą szereg funkcji matematycznych, np. liczbę Pi lub zaokrąglanie
import matplotlib.pyplot as plt  # importuje bibliotekę matplotlib, która jest używana do tworzenia wykresów.
import numpy as np  # importuje bibliotekę numpy, która jest używana do obliczeń numerycznych.
from scipy import \
    signal  # importuje moduł signal z biblioteki scipy, który zawiera funkcje do tworzenia różnych sygnałów.
# importuje funkcje do przekształceń Fouriera z biblioteki scipy.fftpack.
from scipy.fftpack import fft, ifft, fftfreq

# "fft" (Fast Fourier Transform) to funkcja, która pozwala na przekształcenie sygnału z domeny czasowej na domenę
# częstotliwościową. Wynik działania tej funkcji to tzw. widmo sygnału, czyli reprezentacja sygnału w postaci składowych
# częstotliwościowych.

# "ifft" (Inverse Fast Fourier Transform) to funkcja odwrotna do "fft". Służy do przekształcania sygnału z domeny
# częstotliwościowej na domenę czasową.

# "fftfreq" generuje tablicę z odpowiadającymi częstotliwościami dla danej liczby próbek i okresu próbkowania.

R = float(input('Podaj rezystancję(1000=1kOhm): '))  # pobiera od użytkownika wartość rezystancji
C = float(
    input('Podaj pojemność kondensatora(0.00000051=510nF): '))  # pobiera od użytkownika wartość pojemności kondensatora
f = float(input('Podaj częstotliwość(1000=1kHz): '))  # pobiera od użytkownika wartość częstotliwości
T = float(input('Podaj czas pomiaru(0.009=9ms): '))  # pobiera od użytkownika czas pomiaru


def Integrating_System(f, R, C):  # definicja funkcji, która określa działanie układu całkującego
    """Funkcja układu całkującego."""

    # Wejście:
    # F częstotliwość w Hz. (częstotliwość sygnału wejściowego).
    # R wartość rezystancji w Ohmach
    # C wartość pojemności kondensatora w Faradach.

    omega = 2 * np.pi * f  # oblicza wartość omega (częstotliwość kątową) na podstawie podanej częstotliwości f. Dzięki temu omega
    # jest równe ilości pełnych cykli jakie sygnał wykonuje w ciągu sekundy.

    vout = (1. / (
                1j * R * omega * C + 1.))  # oblicza wartość wyjściową układu całkującego. 1./(1j*R*omega*C+1.) jest to wzór transferu
    # układu całkującego. j jest jednostką urojoną, która jest równa sqrt(-1); R jest rezystancją;
    # omega jest częstotliwością kątową, która jest obliczana wyżej; C jest pojemnością kondensatora
    return (vout)  # zwracanie wartości wyjściowej układu


vout_c = Integrating_System(f, R,
                            C)  # wywołanie funkcji układu całkującego i przypisanie jego wartości do zmiennej vout_c

N = 2 ** 15  # ustawia N na 2 do potęgi 15, co jest liczbą punktów próbkowania, które będą wykorzystane przy transformacie Fouriera.
DT = T / N  # obliczanie odstępu między kolejnymi punktami próbkowania, poprzez dzielenie czasu pomiaru przez liczbę punktów
# próbkowania.

t = np.linspace(0., T,
                N)  # tworzenie wektorów czasów, w których próbkowany jest sygnał. np.linspace(0.,T,N) tworzy wektor o
# długości N z równomiernie rozłożonymi punktami od 0 do T.

y_sq = signal.square(
    2 * np.pi * f * t)  # tworzenie sygnału fali prostokątnej o częstotliwości f i czasie t. Funkcja kwadratowa
# pochodzi z scipy.signal i generuje falę prostokątną o podanej częstotliwości i czasie
# w zakresie od 1V do -1V.

f_fft = fftfreq(N, DT)  # obliczanie wartości częstotliwości dla FFT (Fast Fourier Transform) sygnału. fftfreq zwraca
# częstotliwości FFT, biorąc pod uwagę liczbę próbek N i odstęp między próbkami DT.

y_sq_fft = fft(y_sq)  # wykonanie FFT (Fast Fourier Transform) na sygnale prostokątnym zapisanym w y_sq.
# FFT (Fast Fourier Transform) to narzędzie matematyczne, które jest często używane do przekształcania
# sygnału w dziedzinie czasu na odpowiednią reprezentację w dziedzinie częstotliwości. W tym przypadku
# sygnał prostokątny zapisany w y_sq jest przekształcany na jego reprezentację w dziedzinie częstotliwości
# y_sq_fft za pomocą funkcji fft w linii y_sq_fft = fft(y_sq).

y_sq_fft_out = y_sq_fft * Integrating_System(f_fft, R, C)  # użycie systemu całkującego, zdefiniowanego przez funkcję
# Integrating_System do FFT sygnału prostokątnego. Jest to mnożenie
# FFT sygnału prostokątnego przez wynik funkcji Integrating_System.

y_sq_out = ifft(
    y_sq_fft_out)  # wykonanie IFFT (Invert Fast Fourier Transform) na wyjściu systemu całkującego zastosowanego
# w poprzedniej linii. Konwertuje sygnał w dziedzinie częstotliwości z powrotem na sygnał
# w dziedzinie czasu.

plt.figure(figsize=(10, 5))  # tworzy nowy wykres o rozmiarze 10x5

plt.plot(1000 * t, np.real(y_sq), 'black')  # rysowanie pierwszego sygnału (sygnał prostokątny wejściowy) na wykresie z
# kolorem czarnym i zmienienie skali osi x na ms mnożąc czas t przez 10^3

plt.plot(1000 * t, np.real(y_sq_out),
         'r')  # rysowanie drugiego sygnału (sygnał prostokątny po przejściu układu całkującego)
# na wykresie z kolorem czerwonym

ax = plt.gca()  # pobiera aktualny układ wykresu

ax.set_xlim(0.,
            math.ceil(T * 1000))  # ustawia ograniczenie osi x na równy czasowi pomiaru (zaokrąglonemu zawsze w górę)

plt.grid(True)  # dodaje siatkę na wykresie

plt.title("Transient Analysis")  # dodaje tytuł wykresu

plt.xlabel("Time(ms)", position=(0.95, 1))  # dodaje etykietę osi x i ustala jego pozycję

plt.ylabel("Voltage(V)", position=(1, 0.9))  # dodaje etykietę osi y i ustala jego pozycję

plt.show()  # wyswietla wykres