# SZAU — układ zbiorników, część 1

## Uruchomienie

Otwórz katalog projektu w MATLAB-ie R2025b i uruchom `run_part1` w Command Window. Skrypt wymaga **Control System Toolbox**, ale nie wymaga Simulinka. W katalogu `wyniki_czesc1` zapisuje wykresy PNG i cztery osobne pliki PDF z parami wyjścia `y(k)` i sterowania `u(k)`, a tabelę porównania modeli i wskaźniki IAE wyświetla w Command Window. Zapisuje również trzy tabele CSV i plik `dane_sprawozdanie.mat` z wynikami do przygotowania sprawozdania. Wszystkie parametry i przebiegi testowe znajdują się w pliku `run_part1.m`.

## Sprawozdanie

Gotowy dokument: `sprawozdanie_czesc1.pdf` — skrócona wersja ma 8 stron, zawiera podstawowe wzory, sześć rysunków i krótkie opisy wyników. Każdy scenariusz DMC pokazano osobno, ze sterowaniem `u(k)` bezpośrednio pod wyjściem `y(k)`. Edytowalne źródła znajdują się w katalogu `pust_szablon_sprawozdania_utf8`; dokumentem głównym jest `sprawozdanie_szablon.tex`. Szczegóły odtwarzania rysunków i kompilacji zawiera `README_sprawozdanie.md` w tym katalogu. Sprawozdanie obejmuje wyłącznie część 1 projektu.

## Model i punkt pracy

Z treści zadania: `dV1/dt = F1 + FD - α1 sqrt(h1)`, `dV2/dt = α1 sqrt(h1) - α2 sqrt(h2)`, `V1 = A1 h1`, `V2 = C2 h2^3`, `F1(t)=F1in(t-τ)`. Użyte wartości: `A1=560 cm²`, `C2=0,6`, `α1=29`, `α2=33`, `F1in,0=100 cm³/s`, `FD,0=20 cm³/s`, `τ=200 s`. Model nieliniowy całkowany jest funkcją MATLAB-a `ode45` w stanach objętościowych. Przepływy są stałe na kolejnych odcinkach próbkowania o długości 20 s.

Przy `q0=F1in,0+FD,0=120 cm³/s` warunek równowagi daje `α1 sqrt(h1,0)=α2 sqrt(h2,0)=q0`, stąd `h1,0=17,122473 cm` i `h2,0=13,223140 cm`. Pojemność różniczkowa drugiego zbiornika wynosi `B2=3 C2 h2,0²=314,732600 cm²`.

## Model liniowy

Zmienne odchyłkowe: `x=[Δh1; Δh2]`, `Δu=F1in-F1in,0`, `Δd=FD-FD,0`, `Δy=Δh2`. W punkcie pracy `k1=α1/(2 sqrt(h1,0))=3,504167 cm²/s`, `k2=α2/(2 sqrt(h2,0))=4,537500 cm²/s`. Równania:

```text
dx/dt = Ac x + Bc [Δu(t-200); Δd(t)]
Δy = [0 1] x
Ac = [-k1/A1       0       ]
     [ k1/B2      -k2/B2   ]
Bc = [1/A1  1/A1]
     [0     0   ]
```

Transmitancje od sterowania i zakłócenia do wyjścia mają tę samą część właściwą: `G0(s)=k1/[(A1 s+k1)(B2 s+k2)]`; dla sterowania `Gu(s)=G0(s) exp(-200s)`, a dla zakłócenia `Gd(s)=G0(s)`. Wzmocnienie statyczne wynosi `1/k2=0,220385675 cm/(cm³/s)`. Stałe czasowe modelu lokalnego to około 159,8 s i 69,4 s. Skrypt tworzy model przestrzeni stanów funkcją `ss`, dyskretyzuje go przez `c2d` (`Ts=20 s`, opóźnienie 10 próbek), a odpowiedzi na skoki jednostkowe otrzymuje przez `step`. Pozostałe amplitudy uzyskuje przez przemnożenie tych odpowiedzi.

## Porównanie modeli

Skoki sterowania `−50, −25, −10, +10, +25, +50 cm³/s` obejmują wymagane ±50% wartości w punkcie pracy. Analogicznie sprawdzono zakłócenie `FD` w zakresie ±50%, czyli od −10 do +10 cm³/s. Wszystkie próby zaczynają się ze stanu równowagi. Model liniowy i nieliniowy dają identyczny stan początkowy i to samo opóźnienie reakcji na `F1in`; blisko punktu pracy odpowiedzi są zbliżone. Wraz z wielkością skoku różnią się poziom końcowy i kształt przejścia. Dodatnie i ujemne skoki mają różną dynamikę w modelu nieliniowym, gdyż zarówno przepływy pierwiastkowe, jak i pojemność `3 C2 h2²` zależą od poziomu cieczy.

| Skok F1in [cm³/s] | h2 nieliniowy w stanie ustalonym [cm] | h2 liniowy [cm] | Różnica [cm] |
|---:|---:|---:|---:|
| −50 | 4,4995 | 2,2039 | 2,2957 |
| −25 | 8,2874 | 7,7135 | 0,5739 |
| −10 | 11,1111 | 11,0193 | 0,0918 |
| +10 | 15,5188 | 15,4270 | 0,0918 |
| +25 | 19,3067 | 18,7328 | 0,5739 |
| +50 | 26,5381 | 24,2424 | 2,2957 |

Dokładne wartości równowagi po skoku wynikają z `h2,∞=((120+ΔF1in)/33)²`; liniowe z `h2,lin,∞=h2,0+ΔF1in/k2`. Różnica statyczna rośnie kwadratowo z amplitudą skoku. Przy ±50% błąd wynosi już ok. 17% wysokości w punkcie pracy, więc model lokalny nie wystarcza do dokładnego przewidywania pracy w całym badanym zakresie.

## Regulator predykcyjny DMC

Opóźnienie sterowania (200 s) jest dłuższe od szybszej stałej czasowej obiektu. Regulator oparty tylko na bieżącym uchybie reagowałby bez informacji o skutkach poleceń będących jeszcze w torze opóźnienia. DMC przewiduje odpowiedź z uwzględnieniem tej historii. Model pozwala też łatwo rozszerzyć projekt o ograniczenia przepływu w dalszych częściach zadania.

Wartość zadana poziomu to `y_zad(k)`, a wyjście to `y(k)=h2(k)`. Uchyb regulacji `e(k)=y_zad(k)−y(k)` oznacza wartość zadaną minus zmierzony poziom. `Δy_zad` oznacza zmianę zadania względem punktu pracy.

Z liniowej odpowiedzi skokowej `s(i)` zbudowano macierz dynamiczną `M(i,j)=s(i-j+1)` dla `i≥j`, inaczej 0. Swobodna predykcja wynosi `Y0(i)=y(k)+Σ[j=1..D−1] [s(i+j)−s(j)] Δu(k−j)`. Dla wektora wartości zadanych `Y_zad` i funkcji celu `J=||Y_zad−Y0−MΔU||²+λ||ΔU||²` rozwiązanie analityczne to `ΔU=(MᵀM+λI)⁻¹Mᵀ(Y_zad−Y0)`. Stosowany jest tylko pierwszy przyrost sterowania, po czym obliczenia są powtarzane. Parametry: `Ts=20 s`, `D=250`, `N=50`, `Nu=10`, `λ=0,05`. Ta wersja z części 1 nie nakłada ograniczeń na sterowanie.

Testy zamkniętej pętli zaczynają się w tym samym punkcie pracy. Skok wartości zadanej następuje w chwili 400 s, a skok zakłócenia `FD` w chwili 3800 s. Każdy scenariusz uruchamiany jest osobno na obiekcie liniowym i nieliniowym z tym samym regulatorem DMC. IAE dla zadania liczone jest na przedziale 400–3800 s, a dla zakłócenia na przedziale 3800–8000 s. Ostatnia miara oznacza największą różnicę między odpowiedziami obu zamkniętych pętli w całym teście.

| Skok zadania [cm] | Skok FD [cm³/s] | IAE zadanie [cm·s] | IAE zakłócenie [cm·s] | Maks. różnica modeli [cm] |
|---:|---:|---:|---:|---:|
| +2 | +10 | 760 | 873 | 0,44 |
| +4 | +20 | 1714 | 2230 | 1,58 |
| +8 | +30 | 4560 | 4923 | 3,37 |
| −8 | −10 | 3686 | 757 | 2,68 |

Wraz z większymi skokami narasta rozbieżność między modelem liniowym i nieliniowym w zamkniętej pętli. Szczególnie widać to przy zadaniu `h2,0+8 cm`: model nieliniowy ma większe oscylacje. Przy skoku ujemnym `−8 cm` pojawia się oscylacyjna odpowiedź, której model liniowy prawie nie pokazuje. Wszystkie przetestowane wartości `F1in` pozostały dodatnie; najmniejsza wyniosła około 31,2 cm³/s. Regulator pozostaje wersją analityczną bez ograniczeń. Skoki zadania i zakłócenia rosną jednocześnie pomiędzy scenariuszami dodatnimi, dlatego samą różnicę po zakłóceniu należy interpretować jako wynik obu zmian punktu pracy i wielkości `FD`.

Próba `+4 cm, +20 cm³/s` jest zgodna z poprzednim testem. Daje IAE ok. 1714 i 2230 cm·s, a `F1in` pozostaje między około 90 a 134,4 cm³/s. Stanowi punkt odniesienia dla regulatorów rozmytych z części 2.

Wykresy: `wyniki_czesc1/skoki_F1in.png`, `wyniki_czesc1/skoki_FD.png`, `wyniki_czesc1/regulacja_DMC.png` (pierwotny test) i `wyniki_czesc1/regulacja_DMC_porownanie.png` (cztery scenariusze). W MATLAB Workspace pozostają również `P.sys`, `P.sysd` oraz tabele `stepSummary` i `dmcSummary`.
