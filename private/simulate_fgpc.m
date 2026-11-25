function [y,u,decisionTime] = simulate_fgpc(P,model,ctrl,reference,disturbance,isNPL)
% Numeryczny rozmyty regulator predykcyjny FGPC-SL (isNPL=false)
% oraz FGPC-NPL (isNPL=true) z ograniczeniami sterowania u oraz du.
    n = numel(reference)-1;
    y = zeros(n+1,1);
    y(1) = P.h0(2);
    u = zeros(n,1);
    decisionTime = zeros(n,1);
    h = P.h0;

    N = ctrl.N;
    Nu = ctrl.Nu;
    delay = P.delay;
    L = tril(ones(Nu));
    Aineq = [L;-L];
    lb = -ctrl.duMax*ones(Nu,1);
    ub =  ctrl.duMax*ones(Nu,1);
    onesNu = ones(Nu,1);
    lambdaEye = ctrl.lambda*eye(Nu);
    qpOpts = optimoptions('quadprog','Algorithm','active-set','MaxIterations',1000,'Display','off');
    lsqOpts = optimoptions('lsqlin','Display','off');

    for k = 1:n
        if k==1, prevU = P.u0; else, prevU = u(k-1); end
        started = tic;

        % Quasi-liniowy model lokalny w biezacym punkcie pracy y(k).
        [b,c,a0] = local_arx(model,y(k));

        % Odpowiedz skokowa biezacego modelu liniowego -> macierz dynamiczna M_k.
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
            [b0,c0,a00] = local_arx(model,ykm1);
            yMk = b0(1)*ykm1 + b0(2)*ykm2 + c0(1)*u_km1_d + c0(2)*u_km2_d + a00;
        else
            yMk = b(1)*ykm1 + b(2)*ykm2 + c(1)*u_km1_d + c(2)*u_km2_d + a0;
        end
        dk = y(k) - yMk;

        free = zeros(N,1);
        y1 = y(k);
        y2 = ykm1;
        for p = 1:N
            if isNPL
                [bp,cp,a0p] = local_arx(model,y1);
            else
                bp = b; cp = c; a0p = a0;
            end
            u1 = past_or_hold(u,k+p-1-delay,k-1,prevU,P.u0);
            u2 = past_or_hold(u,k+p-2-delay,k-1,prevU,P.u0);
            yNext = bp(1)*y1 + bp(2)*y2 + cp(1)*u1 + cp(2)*u2 + a0p + dk;
            free(p) = yNext;
            y2 = y1;
            y1 = yNext;
        end

        % Zadanie optymalizacji liniowo-kwadratowej z ograniczeniami.
        bineq = [(ctrl.uMax - prevU)*onesNu; (prevU - ctrl.uMin)*onesNu];
        if ctrl.lambda==0
            [moves,~,~,exitflag] = lsqlin(M,reference(k)-free,Aineq,bineq,[],[],lb,ub,[],lsqOpts);
        else
            H = 2*(M'*M + lambdaEye);
            H = 0.5*(H + H');
            f = -2*M'*(reference(k) - free);
            [moves,~,exitflag] = quadprog(H,f,Aineq,bineq,[],[],lb,ub,zeros(Nu,1),qpOpts);
        end
        if exitflag<=0
            error('Optymalizator nie rozwiazal probki k=%d (exitflag=%d).',k,exitflag);
        end
        if any(Aineq*moves>bineq+1e-8) || any(moves<lb-1e-8 | moves>ub+1e-8)
            error('Plan sterowania narusza ograniczenia w probce k=%d.',k);
        end
        decisionTime(k) = toc(started);

        u(k) = prevU + moves(1);
        if k>delay, actualU = u(k-delay); else, actualU = P.u0; end
        h = plant_step(P,h,actualU,disturbance(k));
        y(k+1) = h(2);
    end
end

function [b,c,a0] = local_arx(model,h2)
    w = ts_memberships(h2,model.centers);
    b = model.b*w.';
    c = model.c*w.';
    a0 = model.a0*w.';
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
