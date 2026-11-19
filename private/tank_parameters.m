function P = tank_parameters()
    P.A1 = 560;       % pole przekroju zbiornika 1 [cm^2]
    P.C2 = 0.6;       % V2 = C2*h2^3
    P.a1 = 29;        % wspolczynnik odplywu ze zbiornika 1
    P.a2 = 33;        % wspolczynnik odplywu ze zbiornika 2
    P.u0 = 100;       % F1in w punkcie pracy [cm^3/s]
    P.d0 = 20;        % FD w punkcie pracy [cm^3/s]
    P.tau = 200;      % opoznienie sterowania [s]
    P.Ts = 20;        % okres probkowania [s]
    P.delay = P.tau/P.Ts;
    P.q0 = P.u0+P.d0;
    P.h0 = [(P.q0/P.a1)^2; (P.q0/P.a2)^2];
end
