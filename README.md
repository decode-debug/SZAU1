# SZAU — sterowanie układem dwóch zbiorników

Wzory, wykresy i porównanie regulatorów zawiera [sprawozdanie z części 1–4](sprawozdanie_czesci_1_2_3_4.pdf).

## Symulacje

W MATLAB-ie R2025b otwórz katalog projektu i wpisz **`run_all`**, aby wykonać wszystkie cztery części. Możesz też uruchomić wybraną część osobno. Pełny zestaw wymaga Control System Toolbox, Fuzzy Logic Toolbox i Optimization Toolbox.

| Polecenie w MATLAB-ie | Co oblicza |
|---|---|
| `run_part1` | Odpowiedzi modelu liniowego i nieliniowego oraz regulację klasycznym DMC. |
| `run_part2` | Modele Takagiego–Sugeno z 2–5 regułami i regulację analitycznym PDC przed i po strojeniu. |
| `run_part3` | Regulację numerycznym SL i porównanie z PDC. Wczytuje wyniki części 2. |
| `run_part4` | Regulację numerycznym NPL i porównanie z DMC, PDC oraz SL. Wczytuje wyniki części 3. |

Wspólne parametry zbiorników są w [tank_parameters.m](private/tank_parameters.m), a pętla regulacji SL i NPL w [simulate_fdmc.m](private/simulate_fdmc.m). MATLAB sam znajduje te funkcje w katalogu `private`.

Każdy skrypt wykonuje wszystkie swoje próby i zapisuje wyniki w odpowiednim katalogu od `wyniki_czesc1` do `wyniki_czesc4`. Cztery główne próby regulacji mają pary skoków (wartość zadana w cm, zakłócenie w cm³/s) równe (+2, +10), (+4, +20), (+8, +30) i (−8, −10).

Po zmianie regulatora możesz uruchomić `audit_part3` dla SL lub `audit_part4` dla NPL. Oba testy odtwarzają ruchy sterowania z zapisanych wyników i sprawdzają ograniczenia. Audyt SL sprawdza też m.in. predykcję i wartości IAE.

## Sprawozdanie

Za wygląd i zapis wykresów odpowiada klasa [ResultPlots](@ResultPlots/ResultPlots.m). Jej metody w katalogu `@ResultPlots` przyjmują przebiegi obliczone w `run_part1`–`run_part4`.

W katalogu projektu uruchom `py zbuduj_sprawozdanie.py`. Skrypt korzysta z Pythona 3 i `pdflatex`. Gotowy raport zapisuje jako `sprawozdanie_czesci_1_2_3_4.pdf`.

Po zmianie wyników symulacji odśwież materiały w katalogu `pust_szablon_sprawozdania_utf8` przed kompilacją.

- **Część 1.** Uruchom w tym katalogu `py przygotuj_rysunki.py` (NumPy, SciPy, Matplotlib). Skrypt przygotuje rysunki i tabele.
- **Części 2 i 3.** Skopiuj używane wykresy PNG z `wyniki_czesc2` i `wyniki_czesc3` odpowiednio do `rysunki/czesc2` i `rysunki/czesc3`.
- **Część 4.** LaTeX pobiera wykresy bezpośrednio z `wyniki_czesc4`.

Wartości wpisane ręcznie w rozdziałach i wnioski z prób zaktualizuj w plikach `.tex`.
