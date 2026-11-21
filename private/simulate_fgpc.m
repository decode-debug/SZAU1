function [y,u,decisionTime] = simulate_fgpc(P,model,ctrl,Nu,reference,disturbance,options,isNPL)
% Numeryczny rozmyty regulator predykcyjny FGPC-SL (isNPL=false, wzory 4.103-4.105, 4.107)
% oraz FGPC-NPL (isNPL=true, wzory 4.103-4.105, 4.118-4.125) ze skryptu SZAU.
% Wykorzystuje dyskretny model Takagiego-Sugeno (rownania roznicowe w nastepnikach)
% i optymalizacje kwadratowa quadprog z ograniczeniami sterowania u oraz du.
    n = numel(reference)-1;
    y = zeros(n+1,1);
    y(1) = P.h0(2);
    u = zeros(n,1);
    decisionTime = zeros(n,1);
    h = P.h0;

    N = ctrl.N;
    delay = P.delay;
    L = tril(ones(Nu));
    Aineq = [L;-L];
    lb = -ctrl.duMax*ones(Nu,1);
    ub =  ctrl.duMax*ones(Nu,1);
    onesNu = ones(Nu,1);
    lambdaEye = ctrl.lambda*eye(Nu);

    for k = 1:n
        if k==1, prevU = P.u0; else, prevU = u(k-1); end
        started = tic;

        % (4.107): quasi-liniowy model lokalny w biezacym punkcie pracy y(k).
        w = ts_memberships(y(k),model.centers);
        b = model.b*w.';
        c = model.c*w.';
        a0 = model.a0*w.';

        % Odpowiedz skokowa biezacego modelu liniowego -> macierz dynamiczna M_k (4.101).
        s = zeros(N,1);
        for p = 1:N
            y1 = 0; y2 = 0;
            if p>1, y1 = s(p-1); end
            if p>2, y2 = s(p-2); end
            s(p) = b(1)*y1 + b(2)*y2 + c(1)*(p>delay) + c(2)*(p>delay+1);
        end
        M = zeros(N,Nu);
        for p = 1:N
            mv = 1:min(p,Nu);
            M(p,mv) = s(p-mv+1).';
        end

        % Trajektoria swobodna (przyszle przyrosty sterowania rowne 0, tj. u = prevU).
        ykm1 = y(max(1,k-1));
        ykm2 = y(max(1,k-2));
        u_km1_d = past_or_hold(u,k-1-delay,k-1,prevU,P.u0);
        u_km2_d = past_or_hold(u,k-2-delay,k-1,prevU,P.u0);

        if isNPL
            % (4.118)-(4.119): wyjscie modelu rozmytego w chwili k i skladnik d_k.
            w0 = ts_memberships(ykm1,model.centers);
            b0 = model.b*w0.';
            c0 = model.c*w0.';
            a00 = model.a0*w0.';
            yMk = b0(1)*ykm1 + b0(2)*ykm2 + c0(1)*u_km1_d + c0(2)*u_km2_d + a00;
            dk = y(k) - yMk;

            % (4.120)-(4.125): iteracyjna nieliniowa predykcja odpowiedzi swobodnej.
            free = zeros(N,1);
            y1 = y(k);
            y2 = ykm1;
            for p = 1:N
                wp = ts_memberships(y1,model.centers);
                bp = model.b*wp.';
                cp = model.c*wp.';
                a0p = model.a0*wp.';
                u1 = past_or_hold(u,k+p-1-delay,k-1,prevU,P.u0);
                u2 = past_or_hold(u,k+p-2-delay,k-1,prevU,P.u0);
                yNext = bp(1)*y1 + bp(2)*y2 + cp(1)*u1 + cp(2)*u2 + a0p + dk;
                free(p) = yNext;
                y2 = y1;
                y1 = yNext;
            end
        else
            % FGPC-SL (s. 163): odpowiedz swobodna z zamrozonym modelem liniowym (4.107).
            yMk = b(1)*ykm1 + b(2)*ykm2 + c(1)*u_km1_d + c(2)*u_km2_d + a0;
            dk = y(k) - yMk;
            free = zeros(N,1);
            y1 = y(k);
            y2 = ykm1;
            for p = 1:N
                u1 = past_or_hold(u,k+p-1-delay,k-1,prevU,P.u0);
                u2 = past_or_hold(u,k+p-2-delay,k-1,prevU,P.u0);
                yNext = b(1)*y1 + b(2)*y2 + c(1)*u1 + c(2)*u2 + a0 + dk;
                free(p) = yNext;
                y2 = y1;
                y1 = yNext;
            end
        end

        % (4.103)-(4.105): zadanie optymalizacji liniowo-kwadratowej z ograniczeniami.
        H = 2*(M'*M + lambdaEye);
        H = 0.5*(H + H');
        f = -2*M'*(reference(k) - free);
        bineq = [(ctrl.uMax - prevU)*onesNu; (prevU - ctrl.uMin)*onesNu];
        [moves,~,exitflag] = quadprog(H,f,Aineq,bineq,[],[],lb,ub,zeros(Nu,1),options);
        if exitflag<=0
            error('quadprog nie rozwiazal probki k=%d (exitflag=%d).',k,exitflag);
        end
        if any(Aineq*moves>bineq+1e-8) || any(moves<lb-1e-8 | moves>ub+1e-8)
            error('Plan sterowania narusza ograniczenia w probce k=%d.',k);
        end
        decisionTime(k) = toc(started);

        u(k) = prevU + moves(1);
        if k>delay, actualU = u(k-delay); else, actualU = P.u0; end
        h = rk4_step(@(state) tank_rhs(P,state,actualU,disturbance(k)),h,P.Ts,5);
        h = max(h,0.05);
        y(k+1) = h(2);
    end
end

function val = past_or_hold(u,idx,lastKnown,prevU,u0)
    if idx<1
        val = u0;
    elseif idx<=lastKnown
        val = u(idx);
    else
        val = prevU;
    end
end

function dh = tank_rhs(P,h,u,d)
    f2 = P.a1*sqrt(max(h(1),0));
    f3 = P.a2*sqrt(max(h(2),0));
    cap2 = 3*P.C2*max(h(2),0.05)^2;
    dh = [(u+d-f2)/P.A1; (f2-f3)/cap2];
end
