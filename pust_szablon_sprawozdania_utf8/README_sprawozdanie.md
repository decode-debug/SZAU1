# Sprawozdanie — części 1–3

Dokument główny to `sprawozdanie_szablon.tex`. Opisy kolejnych etapów znajdują się w plikach `czesc2.tex` i `czesc3.tex`. Gotowy dokument łączny: `../sprawozdanie_czesci_1_2_3.pdf`. Starsze pliki PDF z części 1 i 1–2 pozostają w katalogu nadrzędnym.

## Odtworzenie wyników

1. W katalogu nadrzędnym uruchom w MATLAB-ie R2025b kolejno `run_part1`, `run_part2` i `run_part3`. Poprawność wyników części 3 sprawdza `audit_part3`. Do odtworzenia części 2 i 3 potrzebne są także pliki `ts_memberships.m`, `local_step_responses.m`, `dmc_matrices.m` oraz Control System Toolbox, Fuzzy Logic Toolbox i Optimization Toolbox (`quadprog`).
2. Dla rysunków części 1 uruchom tutaj `python przygotuj_rysunki.py` (NumPy, SciPy, Matplotlib).
3. Skopiuj do `rysunki/czesc2/` pliki PNG z `../wyniki_czesc2/`, a do `rysunki/czesc3/` trzy pliki PNG z `../wyniki_czesc3/`.
4. W katalogu sprawozdania uruchom dwukrotnie `pdflatex -interaction=nonstopmode -halt-on-error sprawozdanie_szablon.tex`.

Źródłem liczb i wykresów są skrypty MATLAB oraz pliki CSV i MAT w katalogach `wyniki_czesc1`, `wyniki_czesc2` i `wyniki_czesc3`. Tabelę błędów dla amplitud 10%, 25% i 50% generuje `run_part2`. Część 3 porównuje analityczny PDC (przed i po strojeniu) z trzema wariantami FDMC-SL: o nastawach PDC, przed strojeniem (`N = 50, Nu = 5, λ = 0,2`) oraz po strojeniu (`N = 20, Nu = 3, λ = 0,005, Tf = 160 s`), wraz z weryfikacją na próbach walidacyjnych.
