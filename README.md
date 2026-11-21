# SZAU — sterowanie układem dwóch zbiorników

Wzory, wykresy i porównanie regulatorów zawiera [sprawozdanie z części 1–4](sprawozdanie_czesci_1_2_3_4.pdf).

## Symulacje

W MATLAB-ie R2025b otwórz katalog projektu i wpisz **`run_all`**, aby wykonać wszystkie cztery części. Możesz też uruchomić wybraną część osobno. Pełny zestaw wymaga Control System Toolbox i Optimization Toolbox (`quadprog`).

| Polecenie w MATLAB-ie | Co oblicza |
|---|---|
| `run_part1` | Odpowiedzi modelu liniowego i nieliniowego oraz regulację klasycznym DMC. |
| `run_part2` | Dyskretne modele Takagiego–Sugeno (2–5 reguł, wzory 4.38 i 4.41) oraz analityczny rozmyty DMC (PDC, wzory 4.72–4.75) przed i po strojeniu. |
| `run_part3` | Numeryczny rozmyty regulator predykcyjny FGPC–SL (wzory 4.103–4.105, 4.107), strojenie parametru $\lambda$, walidację oraz porównanie z PDC. Wczytuje wyniki części 2. |
| `run_part4` | Numeryczny rozmyty regulator predykcyjny FGPC–NPL (wzory 4.103–4.105, 4.118–4.125), strojenie parametru $\lambda$, walidację oraz zbiorcze porównanie z DMC, PDC i FGPC–SL. Wczytuje wyniki części 3. |

Wspólne parametry zbiorników są w [tank_parameters.m](private/tank_parameters.m), budowa dyskretnego modelu rozmytego TS w [build_arx_model.m](private/build_arx_model.m), funkcje przynależności w [ts_memberships.m](private/ts_memberships.m), a pętla regulacji FGPC–SL i FGPC–NPL w [simulate_fgpc.m](private/simulate_fgpc.m). MATLAB sam znajduje te funkcje w katalogu `private`.

Każdy skrypt wykonuje wszystkie swoje próby i zapisuje wyniki w odpowiednim katalogu od `wyniki_czesc1` do `wyniki_czesc4`. Cztery główne próby regulacji mają wspólną oś czasu (0–8000 s, skok zadanej w 400 s, skok $F_D$ w 3800 s) i pary skoków (wartość zadana w cm, zakłócenie w cm³/s) równe (+2, +10), (+4, +20), (+8, +30) i (−8, −10).

Po zmianie regulatora możesz uruchomić `audit_part3` dla FGPC–SL lub `audit_part4` dla FGPC–NPL. Oba testy niezależnie odtwarzają ruchy sterowania przez `lsqlin`, sprawdzają ograniczenia, dokładność dyskretyzacji ZOH oraz równoważność symulacji obiektu w zmiennych wysokości i objętości.

## Sprawozdanie

Za wygląd i zapis wykresów odpowiada klasa [ResultPlots](@ResultPlots/ResultPlots.m). Jej metody w katalogu `@ResultPlots` przyjmują przebiegi obliczone w `run_part1`–`run_part4`.

W katalogu projektu uruchom `py zbuduj_sprawozdanie.py`. Skrypt korzysta z Pythona 3 i `pdflatex`. Gotowy raport zapisuje jako `sprawozdanie_czesci_1_2_3_4.pdf`.

Po zmianie wyników symulacji odśwież materiały w katalogu `pust_szablon_sprawozdania_utf8` przed kompilacją.

- **Część 1.** Uruchom w tym katalogu `py przygotuj_rysunki.py` (NumPy, SciPy, Matplotlib). Skrypt przygotuje rysunki i tabele.
- **Części 2 i 3.** Skopiuj używane wykresy PNG z `wyniki_czesc2` i `wyniki_czesc3` odpowiednio do `rysunki/czesc2` i `rysunki/czesc3`.
- **Część 4.** LaTeX pobiera wykresy bezpośrednio z `wyniki_czesc4`.

Wartości wpisane ręcznie w rozdziałach i wnioski z prób zaktualizuj w plikach `.tex`.
