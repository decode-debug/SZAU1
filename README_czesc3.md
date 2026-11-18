# SZAU — układ zbiorników, część 3

## Uruchomienie

Z katalogu projektu uruchom w MATLAB-ie R2025b kolejno `run_part2`, `run_part3` i `audit_part3`. Obliczenia używają Control System Toolbox, Fuzzy Logic Toolbox i Optimization Toolbox. Część 3 wczytuje czteroregułowy model oraz przebiegi PDC z `wyniki_czesc2/dane_czesc2.mat`.

## FDMC-SL według materiałów wykładowych

Podstawą są wzory (4.96)–(4.105) na s. 161–162 skryptu „Regulatory rozmyte”. Funkcja `local_step_responses.m` oblicza przed symulacją znormalizowaną odpowiedź poziomu na jednostkowy skok dopływu dla każdego modelu lokalnego. Kolumna `S(:,i)` zawiera próbki od 1 do D, łącznie z opóźnieniem 200 s. PDC i SL korzystają z tych samych lokalnych odpowiedzi.

W chwili k funkcja `ts_memberships.m` wyznacza wagi na podstawie pomiaru h₂ i węzłów modelu po strojeniu. FDMC-SL oblicza `s = S*w` dla kolumnowego wektora znormalizowanych wag. Wagi są stałe podczas tej predykcji. Funkcja `dmc_matrices.m` buduje z tych samych rzędnych macierz dynamiczną i macierz wpływu przeszłych ruchów. Po D próbkach przyjmuje stałą wartość ostatniej rzędnej.

`quadprog` minimalizuje sumę kwadratów błędów prognozy i karę λ za planowane przyrosty. Dla każdej planowanej próbki obowiązuje 0 ≤ u ≤ 200 cm³/s i |Δu| ≤ 30 cm³/s. W dostrojonym wariancie SL wartość zadana podawana do funkcji celu przechodzi przez przyczynowy filtr pierwszego rzędu o stałej czasowej Tf = 160 s (`alpha = exp(-Ts/Tf)`), który korzysta wyłącznie z bieżącej wartości zadanej. Regulator wykonuje pierwszy ruch i wyznacza nową predykcję po kolejnym pomiarze.

## Warunki porównania

Wszystkie regulatory sterują tym samym nieliniowym obiektem przy D = 250 i okresie próbkowania 20 s. Skok zadania występuje w 400 s, skok zakłócenia w 3200 s, a symulacja kończy się w 6400 s.

| Regulator | Horyzont predykcji N | Horyzont sterowania Nu | Kara λ | Filtr Tf [s] | Wagi |
|---|---:|---:|---:|---:|---|
| Analityczny PDC | 50 | 1 | 0,05 | 0 | Jak w modelu TS |
| Analityczny PDC po strojeniu | 50 | 1 | 0,05 | 0 | Węzły regulatora i γ = 1,5 z części 2 |
| FDMC-SL o nastawach PDC | 50 | 1 | 0,05 | 0 | Jak w modelu TS |
| FDMC-SL przed strojeniem | 50 | 5 | 0,2 | 0 | Jak w modelu TS |
| FDMC-SL po strojeniu | 20 | 3 | 0,005 | 160 | Jak w modelu TS |

Pierwszy PDC i pierwszy SL mają identyczne nastawy, ograniczenia i funkcje przynależności. PDC uśrednia dopuszczalne ruchy lokalnych regulatorów. SL uśrednia lokalne odpowiedzi skokowe i rozwiązuje jedno wspólne zadanie kwadratowe.

## Wyniki

| Skok zadania / zakłócenia | PDC | PDC strojony | SL o nastawach PDC | SL przed strojeniem | SL po strojeniu |
|---|---:|---:|---:|---:|---:|
| +2 cm / +10 cm³/s | 2082 | 1861 | 2085 | 1919 | 1676 |
| +4 cm / +20 cm³/s | 4568 | 4250 | 4566 | 4159 | 3572 |
| +8 cm / +30 cm³/s | 9285 | 8718 | 9181 | 8325 | 6867 |
| −8 cm / −10 cm³/s | 3601 | 3478 | 3656 | 3516 | 3685 |
| **Łączny IAE [cm·s]** | **19 535** | **18 309** | **19 488** | **17 919** | **15 800** |

Przy identycznych nastawach SL ma łączny IAE niższy od PDC o 0,24%. Wariant SL przed strojeniem (`N = 50, Nu = 5, λ = 0,2`) obniża łączny IAE o 2,13% względem dostrojonego PDC, ale przed osiągnięciem nowej wartości zadanej wykazuje wyraźne cofnięcie poziomu (do 0,575 cm w próbie +8 cm) oraz niewielkie przeregulowanie (do 0,073 cm).

Skorygowanie nastaw do `N = 20, Nu = 3, λ = 0,005` wraz z przyczynowym filtrem wartości zadanej `Tf = 160 s` obniża łączny IAE w czterech próbach projektowych do 15 800 cm·s (o 11,82% mniej niż SL przed strojeniem i o 13,70% mniej niż dostrojony PDC). IAE w fazie tłumienia zakłócenia spada z 8260 do 5471 cm·s (o 33,76%), a maksymalne odchylenie po skoku FD +30 cm³/s maleje z 4,90 do 3,95 cm. Przeregulowanie znika we wszystkich próbach (0,000 cm), a maksymalne cofnięcie na dojściu do zadania maleje niemal dziewięciokrotnie (z 0,575 do 0,066 cm; w próbie +8 cm do 0,023 cm). Maksymalny pojedynczy przyrost sterowania spada z 15,60 do 11,72 cm³/s przy sumie zmian 310,1 cm³/s.

Na czterech dodatkowych próbach walidacyjnych o innych amplitudach i przeciwnych znakach wymuszeń (`+3/−8`, `−3/+12`, `+6/−10`, `−6/+15`) wariant SL po strojeniu wygrywa w każdym scenariuszu, obniżając łączny IAE z 12 292 cm·s (SL o nastawach PDC) i 11 523 cm·s (SL przed strojeniem) do 11 100 cm·s (o 3,67% mniej niż przed strojeniem), przy zerowym przeregulowaniu i cofnięciu ograniczonym z 0,457 do 0,058 cm.

Pliki `porownanie_wyjsc.png`, `porownanie_sterowan.png` i `strojenie_SL.png` pokazują przebiegi. `porownanie_regulatorow.csv` zawiera wskaźniki pięciu regulatorów, `wskazniki_SL.csv` nastawy i czasy decyzji trzech wariantów SL na próbach projektowych, a `walidacja_SL.csv` wyniki na próbach walidacyjnych.

## Sprawdzenie wyników

`audit_part3.m` porównuje lokalne odpowiedzi z niezależną dyskretyzacją ZOH, predykcję dla dwóch równocześnie aktywnych reguł z równoległą symulacją ich stanów oraz każdy ruch QP we wszystkich wariantach SL (dla próby projektowej i walidacyjnej) z niezależną rekonstrukcją celu i rozwiązaniem `lsqlin`. Sprawdza także wagi, ograniczenia zrealizowanych ruchów, cofnięcie, przeregulowanie, IAE i bilans objętości.
