# SZAU — układ zbiorników, część 2

## Uruchomienie

W MATLAB-ie R2025b otwórz katalog projektu i wykonaj `run_part2`. Wymagane są **Control System Toolbox i Fuzzy Logic Toolbox**. Funkcje przynależności oblicza `trimf`; lokalne modele dynamiczne i analityczne prawo DMC pozostają jawne w kodzie. MPC Toolbox nie jest potrzebny. Wyniki liczbowe pojawiają się w Command Window, a pliki PNG, CSV i MAT w katalogu `wyniki_czesc2`.

W pliku `run_part2.m` sekcje `%%` odpowiadają kolejno: parametrom, modelom TS, ocenie modeli, strojeniu modelu, projektowi DMC, strojeniu regulatora i porównaniu zamkniętych pętli. Wszystkie testy można odtworzyć przez jedno uruchomienie. W Workspace pozostają m.in. `modelSummary`, `controlSummary`, `modelTuned`, `controllerCenters` i `traces`.

## Modele Takagiego–Sugeno

Wybrano poziom `h2` jako zmienną reguł, ponieważ jest mierzonym wyjściem, a pojemność różniczkowa drugiego zbiornika `3 C2 h2²` mocno zmienia się wraz z tym poziomem. Zakres funkcji przynależności to `h2 = 4,4995..26,5381 cm`; przy `FD=20 cm³/s` odpowiada on dopływom `F1in=50..150 cm³/s` (±50% punktu pracy). Podział początkowy jest równomierny, z trójkątnymi funkcjami i ramionami na brzegach. Wagi są nieujemne i sumują się do 1.

W każdym węźle `h2,i` wyznaczono lokalny punkt równowagi `qi=α2 sqrt(h2,i)` i `h1,i=(qi/α1)²`, a następnie macierz Jacobiego:

```text
ki1 = α1/(2 sqrt(h1,i)),  ki2 = α2/(2 sqrt(h2,i))
Bi2 = 3 C2 h2,i²
Ai  = [-ki1/A1      0      ]
      [ ki1/Bi2   -ki2/Bi2 ]
B   = [1/A1; 0]
fi(h,u,d) = Ai (h-hi) + B (u+d-qi)
fTS(h,u,d) = Σ μi(h2) fi(h,u,d)
```

Stanami są `[h1;h2]`, `u` oznacza rzeczywisty dopływ do zbiornika po opóźnieniu, a `d=FD`. Symulacja TS używa tego samego okresu 20 s i opóźnienia 200 s co model nieliniowy. Ponieważ lokalne modele są liniaryzacjami w *różnych* punktach równowagi, ich średnia nie musi dokładnie zachowywać punktu równowagi pośredniego; jest to jedna z przyczyn błędu dla małej liczby reguł.

Ocena obejmuje 12 odpowiedzi od punktu pracy: skoki `F1in = ±10, ±25, ±50 cm³/s` oraz `FD = ±2, ±5, ±10 cm³/s`. RMSE liczono po wszystkich próbkach tych przebiegów. Model liniowy z części 1 stanowi punkt odniesienia:

| Model | RMSE [cm] | Największy błąd [cm] |
|---|---:|---:|
| Liniowy (1 lokalny) | 0,878 | 2,296 |
| TS, 2 reguły | 2,213 | 3,834 |
| TS, 3 reguły | 0,771 | 2,490 |
| TS, 4 reguły | 0,220 | 1,224 |
| TS, 5 reguł | 0,129 | 0,703 |

Model z dwiema regułami jest zbyt zgrubny i na tym zestawie bywa gorszy od modelu liniowego. Wybrano model czteroregułowy: poprawia RMSE ponad trzykrotnie względem modelu trzyregułowego, zachowując o jedną regułę mniej niż najlepszy pięcioregułowy.

Po pogrupowaniu skoków według procentowej zmiany każdego wejścia na grupę przypadają cztery odpowiedzi (oba wejścia, oba znaki). Przy 10% najdokładniejszy jest model liniowy z RMSE 0,060 cm. Przy 50% jego RMSE rośnie do 1,472 cm, a dla TS4 wynosi 0,201 cm i dla TS5 0,063 cm. Zestawienie dla wszystkich trzech amplitud generuje `jakosc_modeli_amplitudy.csv`. Przy stałych dopływach w pierwszych 200 s poziom w TS2 maleje o 1,854 cm, a w TS4 o 0,066 cm. Zjawisko wynika z dryfu modelu między punktami linearyzacji i występuje także w czasie opóźnienia skoku sterowania.

### Strojenie modelu

Dwa wewnętrzne węzły funkcji przynależności modelu czteroregułowego przesunięto na siatce od −3 do +3 cm co 1,5 cm, pozostawiając końce zakresu i **macierze lokalnych modeli** stałe. Początkowe węzły `[4,500; 11,846; 19,192; 26,538]` zmieniły się na `[4,500; 10,346; 19,192; 26,538]` cm. Na 12 skokach użytych do doboru RMSE spadło z 0,220 do 0,185 cm (o około 16%). Korzyścią jest lepsze dopasowanie w badanym zestawie przy tej samej liczbie reguł i tych samych lokalnych linearyzacjach.

Sprawdzono także cztery dodatkowe skoki (`F1in=±40`, `FD=±7 cm³/s`), niewykorzystane do wyboru węzłów. Tam RMSE **wzrosło** z 0,167 do 0,175 cm. Strojenie nie poprawia więc każdego wymuszenia; ten wynik należy uwzględnić przy ocenie przydatności modelu poza zestawem projektowym.

## Analityczny rozmyty DMC z ograniczeniami

Na podstawie wybranego modelu TS wyznaczono dla każdej reguły lokalną odpowiedź skokową i lokalny regulator DMC. Dla regulacji z ograniczeniami przyjęto `N=50`, `D=250`, `Nu=1`, `λ=0,05`, `Ts=20 s`. Jedna optymalizowana zmiana sterowania (`Nu=1`) pozwala policzyć dokładne rozwiązanie analityczne przy ograniczeniach bieżącego ruchu:

```text
y0,i = y(k)·1 + Fi·Δu_past
gi   = [si(1), ..., si(N)]ᵀ
Δui* = giᵀ(r·1 − y0,i) / (giᵀgi + λ)
L    = max(−ΔuMax, uMin − u(k−1))
H    = min( ΔuMax, uMax − u(k−1))
Δui  = min(max(Δui*,L),H)
Δu   = Σ μi(h2(k)) Δui
```

Użyto `0 ≤ F1in ≤ 200 cm³/s` i `|ΔF1in| ≤ 30 cm³/s` na próbkę. Każdy lokalny ruch jest rzutowany na ten sam przedział `[L,H]`, więc ich suma ważona również spełnia oba ograniczenia. Ograniczenia dotyczą bieżącego sterowania i jego zmiany; w tej wersji nie ograniczano przewidywanego wyjścia. Wagi regulatora PDC mają początkowo ten sam kształt co w wybranym modelu TS. Model obiektu w zamkniętej pętli zawsze jest nieliniowy.

Do porównania zbudowano też **liniowy DMC z tymi samymi** `N, D, Nu, λ` i ograniczeniami. Osobno pokazano regulator liniowy z części 1 (`Nu=10`, bez ograniczeń), dlatego jego wyniki należy interpretować z uwzględnieniem innego horyzontu sterowania.

### Testy i strojenie regulatora

Każda próba zaczyna się w punkcie pracy. W chwili 400 s zmienia się wartość zadana, w 3200 s dopływ zakłócający. Wskaźniki IAE liczono odpowiednio na przedziałach 400–3200 s i 3200–6400 s. Regulator rozmyty strojono przez potęgowanie wag i ich ponowną normalizację (`μi^γ/Σμj^γ`) oraz przesunięcie dwóch wewnętrznych węzłów *tylko w regulatorze*. Model obiektu i lokalne odpowiedzi skokowe pozostały te same. Najlepszy wynik na czterech scenariuszach dało `γ=1,5` i węzły `[4,500; 12,346; 21,192; 26,538]` cm, wobec początkowych `[4,500; 10,346; 19,192; 26,538]` cm.

| Skok zadania / FD | DMC część 1 | Liniowy z ograniczeniami | Rozmyty PDC | Rozmyty po strojeniu |
|---|---:|---:|---:|---:|
| +2 cm / +10 cm³/s | 1633 | 1890 | 2082 | 1861 |
| +4 cm / +20 cm³/s | 3937 | 4425 | 4568 | 4250 |
| +8 cm / +30 cm³/s | 9285 | 10148 | 9285 | 8718 |
| −8 cm / −10 cm³/s | 4378 | 4953 | 3601 | 3478 |
| **Suma IAE [cm·s]** | **19232** | **21416** | **19535** | **18309** |

W tabeli dla każdego scenariusza zsumowano IAE po zmianie zadania i zakłócenia. Po strojeniu suma IAE spadła o około 6,3% względem początkowego PDC oraz o około 14,5% względem liniowego regulatora z tymi samymi ograniczeniami. Na większych skokach rozmyty regulator zyskuje najwięcej; dla małego skoku +2 cm DMC z części 1 pozostaje szybszy. Sumarycznie dostrojony regulator rozmyty osiąga o około 4,8% niższy IAE także od DMC z części 1, choć ten drugi nie ma narzuconych ograniczeń. We wszystkich ograniczonych próbach odnotowano `0 ≤ F1in ≤ 200` i `|ΔF1in| ≤ 30`.

Wyniki strojenia regulatora dotyczą scenariuszy użytych do wyboru funkcji przynależności; poprawy dla innych wymuszeń nie należy zakładać bez osobnego testu.

## Pliki wynikowe

- `funkcje_przynaleznosci.png` i `funkcje_przynaleznosci_strojone.png`: funkcje modeli TS.
- `modele_2_5.png` i `porownanie_modeli_TS.png`: odpowiedzi modeli przed i po strojeniu.
- `regulacja_porownanie.png` i `strojenie_regulatora.png`: działanie regulatorów i zmiana wag.
- `jakosc_modeli.csv`, `jakosc_modeli_amplitudy.csv`, `jakosc_regulacji.csv` i `dane_czesc2.mat` zawierają wskaźniki i przebiegi do sprawozdania.
