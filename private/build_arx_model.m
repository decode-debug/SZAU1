function model = build_arx_model(P,centers)
% Dyskretny dynamiczny model rozmyty Takagiego-Sugeno.
% Dla kazdego wezla h2_i wyznacza lokalna linearyzacje ukladu dwoch zbiornikow
% i jej dokladna dyskretyzacje ZOH w postaci rownania roznicowego:
%   y^i(k+1) = b1_i*y(k) + b2_i*y(k-1) + c1_i*q(k) + c2_i*q(k-1) + a0q_i,
% gdzie q(k) = F1in(k-delay) + FD(k). Przy nominalnym FD = d0:
%   y^i(k+1) = b1_i*y(k) + b2_i*y(k-1) + c1_i*u(k-delay) + c2_i*u(k-delay-1) + a0_i.
    model.centers = centers(:)';
    model.eqCenters = centers(:)';
    model.count = numel(centers);
    model.A = zeros(2,2,model.count);
    model.B = [1/P.A1;0];
    model.heq = zeros(2,model.count);
    model.qeq = zeros(1,model.count);
    model.b = zeros(2,model.count);
    model.c = zeros(2,model.count);
    model.a0q = zeros(1,model.count);
    model.a0 = zeros(1,model.count);

    for i = 1:model.count
        h2 = centers(i);
        q = P.a2*sqrt(h2);
        h1 = (q/P.a1)^2;
        k1 = P.a1/(2*sqrt(h1));
        k2 = P.a2/(2*sqrt(h2));
        cap2 = 3*P.C2*h2^2;
        A = [-k1/P.A1, 0; k1/cap2, -k2/cap2];

        model.heq(:,i) = [h1;h2];
        model.qeq(i) = q;
        model.A(:,:,i) = A;

        sysd = c2d(ss(A,model.B,[0 1],0),P.Ts,'zoh');
        [num,den] = tfdata(sysd,'v');
        b1 = -den(2);
        b2 = -den(3);
        c1 = num(2);
        c2 = num(3);
        model.b(:,i) = [b1;b2];
        model.c(:,i) = [c1;c2];
        model.a0q(i) = (1-b1-b2)*h2 - (c1+c2)*q;
        model.a0(i)  = (1-b1-b2)*h2 - (c1+c2)*(q-P.d0);
    end
end
