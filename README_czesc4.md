# Część 4 — FDMC-NPL

W MATLAB-ie uruchom z katalogu projektu kolejno `run_part2`, `run_part3`, `run_part4` i `audit_part4`. Jeśli masz aktualny `wyniki_czesc3/dane_czesc3.mat`, wystarczą dwa ostatnie polecenia. Potrzebne są te same pakiety co w części 3.

`simulate_fdmc.m` realizuje wspólną pętlę SL i NPL. Podanie dodatkowego argumentu z modelem TS włącza predykcję NPL. Lokalny poziom w chwili pomiaru wynika z historii przyrostów dopływu, wzmocnienia ustalonego oraz lokalnego punktu równowagi. Regulator oblicza poprawkę pomiarową z wagami z poprzedniej próbki. Na horyzoncie utrzymuje tę poprawkę, zakłada zerowe przyszłe ruchy i aktualizuje wagi na podstawie ostatniego przewidzianego poziomu. Macierz dynamiczną buduje z bieżącego pomiaru tak jak SL.

W katalogu `wyniki_czesc4` znajdują się miary każdego scenariusza, podsumowanie, przebiegi MAT i wykresy wyjścia oraz sterowania. Wiersze ośmiu regulatorów odnoszą się do tych samych czterech prób. Audyt odtwarza decyzje NPL przez bezpośredni splot odpowiedzi impulsowych z historią sterowania oraz niezależne rozwiązanie zadania najmniejszych kwadratów. Sprawdza też ograniczenia wykonanych ruchów.

Pełny raport to `sprawozdanie_czesci_1_2_3_4.pdf`. Aby go odtworzyć po symulacji, w katalogu `pust_szablon_sprawozdania_utf8` skompiluj dwukrotnie `sprawozdanie_1_2_3_4.tex`.
