from pathlib import Path

root = Path(__file__).resolve().parent.parent
folders = [root/'czesc1_sprawozdanie'/'pust_szablon_sprawozdania_utf8',
           root/'pust_szablon_sprawozdania_utf8']

content = {
'wstep.tex': r'''\chapter{Model obiektu}
Zbadano odpowiedzi układu dwóch zbiorników, wyznaczono model liniowy w punkcie pracy i sprawdzono działanie konwencjonalnego regulatora DMC. Wyjściem jest poziom $h_2$, sterowaniem dopływ $\uvar$, a zakłóceniem dopływ $F_D$.

\section{Model nieliniowy i punkt pracy}
Model opisują bilanse objętości oraz zależności geometryczne:
\begin{align}
\dot V_1&=\uvar(t-\tau)+F_D-\alpha_1\sqrt{h_1},\\
\dot V_2&=\alpha_1\sqrt{h_1}-\alpha_2\sqrt{h_2},\\
V_1&=A_1h_1,\qquad V_2=C_2h_2^3.
\end{align}
Przyjęto $A_1=560\ \mathrm{cm^2}$, $C_2=\num{0.6}$, $\alpha_1=29$, $\alpha_2=33$ oraz $\tau=200\ \mathrm{s}$. W punkcie pracy $\uop=100\ \jednprzep$ i $F_{D,0}=20\ \jednprzep$. Warunek równowagi daje
\begin{equation}
h_{1,0}=\left(\frac{120}{29}\right)^2\approx\num{17.1225}\ \mathrm{cm},\qquad
h_{2,0}=\left(\frac{120}{33}\right)^2\approx\num{13.2231}\ \mathrm{cm}.
\end{equation}
Wszystkie próby rozpoczynano od tego samego punktu równowagi. Zastosowano okres próbkowania $T_s=20\ \mathrm{s}$; opóźnienie sterowania odpowiada 10 próbkom. Zakłócenie $F_D$ działa bez tego opóźnienia.
''',
'wzory.tex': r'''\section{Model liniowy}
Zlinearyzowano równania względem poziomów w punkcie pracy. Dla $x=[\Delta h_1,\Delta h_2]^{\mathsf T}$ otrzymano
\begin{equation}
\dot x=
\begin{bmatrix}-\num{0.00625744}&0\\\num{0.01113379}&-\num{0.01441700}\end{bmatrix}x+
\begin{bmatrix}1/560&1/560\\0&0\end{bmatrix}
\begin{bmatrix}\Delta u(t-200)\\\Delta d(t)\end{bmatrix},\quad
\Delta y=\begin{bmatrix}0&1\end{bmatrix}x,
\end{equation}
gdzie $\Delta u=\uvar-100$, $\Delta d=F_D-20$ oraz $\Delta y=h_2-h_{2,0}$. Transmitancje mają postać
\begin{equation}
G_u(s)=\frac{\num{0.220386}\,e^{-200s}}{(\num{159.81}s+1)(\num{69.36}s+1)},\qquad
G_d(s)=G_u(s)e^{200s}.
\end{equation}
Model liniowy ma dwie inercje i dodatnie wzmocnienie statyczne. Różnica między torami sterowania i zakłócenia polega na opóźnieniu.
''',
'tabele.tex': r'''\chapter{Porównanie modeli}
Porównano odpowiedzi obu modeli na skoki wejść do $\pm50\%$ wartości nominalnych. Każda próba rozpoczyna się od podanego punktu równowagi.
\begin{figure}[H]
\centering
\includegraphics[width=0.90\textwidth]{rysunki/czesc1/skoki_F1in.pdf}
\caption{Skoki sterowania o $\pm10$, $\pm25$ i $\pm50\ \jednprzep$ przy stałym $F_D=20\ \jednprzep$.}
\label{fig:skokiU}
\end{figure}
Blisko punktu pracy odpowiedzi są zbliżone. Dla większych skoków różnią się poziom końcowy i szybkość zmian. Błąd statyczny rośnie jak $(\Delta\uvar)^2/33^2$: dla skoków $\pm10\ \jednprzep$ wynosi około $\num{0.092}\ \mathrm{cm}$, a dla $\pm50\ \jednprzep$ około $\num{2.30}\ \mathrm{cm}$.
\begin{figure}[H]
\centering
\includegraphics[width=0.90\textwidth]{rysunki/czesc1/skoki_FD.pdf}
\caption{Skoki zakłócenia o $\pm2$, $\pm5$ i $\pm10\ \jednprzep$ przy stałym $\uvar=100\ \jednprzep$.}
\label{fig:skokiD}
\end{figure}
Zakłócenie działa bez opóźnienia $200\ \mathrm{s}$. Mniejsze różnice modeli wynikają z mniejszych amplitud bezwzględnych: skok $50\%$ zakłócenia to $10\ \jednprzep$. Nieliniowa dynamika zmienia się z poziomem cieczy, czego stały model liniowy nie odwzorowuje.
''',
'rysunki.tex': r'''\chapter{Regulacja DMC}
\section{Algorytm i testy}
Zaprojektowano analityczny DMC z modelu liniowego. Uwzględnia on wcześniejsze sterowania, co jest istotne przy opóźnieniu $200\ \mathrm{s}$. Dla odpowiedzi skokowej $s_i$ przyjęto $M_{ij}=s_{i-j+1}$ dla $i\geq j$, a poza tym zero. Swobodna predykcja wynosi
\begin{equation}
Y_i^0(k)=y(k)+\sum_{j=1}^{D-1}(s_{i+j}-s_j)\Delta u(k-j).
\end{equation}
Minimalizowano błąd predykcji i zmiany sterowania:
\begin{align}
J&=\|R-Y^0-M\Delta U\|^2+\lambda\|\Delta U\|^2,\\
\Delta U&=(M^{\mathsf T}M+\lambda I)^{-1}M^{\mathsf T}(R-Y^0).
\end{align}
W każdej próbce stosowano pierwszy przyrost $\Delta U$. Nastawy: $D=250$, $N=50$, $N_u=10$, $\lambda=\num{0.05}$, $T_s=20\ \mathrm{s}$. Regulator działa bez ograniczeń; zakłócenie tłumi przez sprzężenie zwrotne.

Ten sam regulator sprawdzono na obu obiektach. W chwili $400\ \mathrm{s}$ zmieniano zadanie o $\Delta r$, a w chwili $3800\ \mathrm{s}$ zakłócenie o $\Delta F_D$. Amplitudy czterech prób podano na wykresach.
\begin{figure}[H]
\centering
\includegraphics[width=0.90\textwidth]{rysunki/czesc1/DMC_wyjscie.pdf}
\caption{Śledzenie zadania i tłumienie zakłócenia dla czterech scenariuszy. Pionowa linia oznacza chwilę zmiany $F_D$.}
\label{fig:dmcY}
\end{figure}
Przy małych skokach odpowiedzi są zbliżone. Dla $\Delta r=\pm8\ \mathrm{cm}$ obiekt nieliniowy wykazuje wyraźniejsze oscylacje. Maksymalna różnica wyjść dla skoku $+8\ \mathrm{cm}$ wynosi około $\num{3.37}\ \mathrm{cm}$.

\clearpage
\section{Sterowanie i jakość regulacji}
Porównano również dopływy wyznaczane przez regulator. Przy zmianie zakłócenia sterowanie koryguje się w przeciwnym kierunku, lecz jego wpływ dociera do obiektu po czasie transportu.
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{rysunki/czesc1/DMC_sterowanie.pdf}
\caption{Sterowanie w tych samych czterech próbach. Linie ciągłe: obiekt nieliniowy; przerywane: liniowy.}
\label{fig:dmcU}
\end{figure}
Jakość regulacji obiektu nieliniowego oceniono wskaźnikiem
\begin{equation}
\mathrm{IAE}=T_s\sum_k|r(k)-y(k)|.
\end{equation}
Dla zadania sumowano próbki od $400$ do chwili poprzedzającej $3800\ \mathrm{s}$, a dla zakłócenia od $3800$ do $8000\ \mathrm{s}$.
\begin{table}[H]
\centering\small
\caption{Wyniki regulacji obiektu nieliniowego i różnica wyjść obu pętli}
\begin{tabular}{rrrrr}
\toprule
$\Delta r$ & $\Delta F_D$ & $\mathrm{IAE}_r$ & $\mathrm{IAE}_d$ & $\max|y_{\mathrm{NL}}-y_{\mathrm{L}}|$\\
{}[cm] & [$\jednprzep$] & [$\mathrm{cm\cdot s}$] & [$\mathrm{cm\cdot s}$] & [cm]\\
\midrule
\tableinput{dane/DMC_wskazniki.tex}
\bottomrule
\end{tabular}
\end{table}
W badanych próbach sterowanie pozostało dodatnie; najmniejsza wartość wynosiła około $\num{31.19}\ \jednprzep$. Przy większych skokach rosną oscylacje i rozbieżności modeli. W dodatnich scenariuszach zwiększano jednocześnie zmianę zadania i zakłócenia, dlatego sam wzrost IAE nie pozwala oddzielić wpływu tych dwóch czynników.

\clearpage
\section{Przykład: skok zadania o 4 cm}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{rysunki/czesc1/DMC_przypadek_4cm.pdf}
\caption{Poziom cieczy, sterowanie i zakłócenie dla $\Delta r=+4\ \mathrm{cm}$ oraz $\Delta F_D=+20\ \jednprzep$.}
\label{fig:dmc4}
\end{figure}
Po zmianie zadania poziom dochodzi w jego pobliże. Dodatni skok zakłócenia powoduje przejściowy wzrost $h_2$, po którym regulator zmniejsza dopływ i przywraca poziom w pobliże zadania. Sterowanie obiektu nieliniowego mieści się w zakresie około $90$--$\num{134.4}\ \jednprzep$.
''',
'uwagi.tex': r'''\section{Wnioski}
\begin{itemize}
\item Model liniowy dobrze przybliża obiekt w pobliżu punktu pracy. Przy dużych zmianach dopływu błędy poziomu i dynamiki stają się wyraźne.
\item DMC uwzględnia opóźnienie i zapewnia śledzenie zadania oraz tłumienie zakłóceń w badanych próbach. Przy zmianach zadania o $\pm8\ \mathrm{cm}$ odpowiedź nieliniowa jest bardziej oscylacyjna.
\item Jeden model lokalny nie opisuje równie dobrze całego zakresu pracy. Wyniki stanowią punkt odniesienia dla regulatorów rozmytych w dalszych częściach projektu.
\end{itemize}
\smallskip
{\small Źródła: treść zadania \texttt{main (50).pdf}, część 1; skrypt \texttt{run\_part1.m}.}
'''
}

for folder in folders:
    if not folder.is_dir():
        continue
    main = folder/'sprawozdanie_szablon.tex'
    text = main.read_text(encoding='utf-8')
    text = text.replace('\\pagenumbering{roman}\n\\tableofcontents\n\\clearpage\n','')
    text = text.replace('\\input{listingi}\n','').replace('\\input{literatura}\n','')
    if '\\raggedbottom' not in text:
        text = text.replace('\\frenchspacing','\\frenchspacing\n\\raggedbottom')
    main.write_text(text,encoding='utf-8')
    for name, body in content.items():
        (folder/name).write_text(body,encoding='utf-8')
print('Skrócono istniejące źródła; zachowano pięć wykresów.')
