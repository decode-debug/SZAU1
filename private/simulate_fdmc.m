function [y,u,decisionTime] = simulate_fdmc(P,localSteps,centers,ctrl,Nu,Tf,reference,disturbance,options,model)
    n = numel(reference)-1;
    y = zeros(n+1,1); y(1) = P.h0(2);
    u = zeros(n,1);
    decisionTime = zeros(n,1);
    h = P.h0;
    pastDu = zeros(ctrl.D-1,1);
    cumulative = tril(ones(Nu));
    Aineq = [cumulative;-cumulative];
    lb = -ctrl.duMax*ones(Nu,1);
    ub = ctrl.duMax*ones(Nu,1);
    target = reference(1);
    alpha = 0;
    if Tf>0, alpha = exp(-P.Ts/Tf); end

    nonlinearFree = nargin==10;
    if nonlinearFree
        for i = 1:model.count
            [~,localFree(:,:,i)] = dmc_matrices(localSteps(:,i),ctrl.N,Nu);
        end
        gain = localSteps(end,:);
        equilibriumU = model.qeq-P.d0;
        historyGain = localSteps(1:end-1,:)-gain;
    end
    for k = 1:n
        if k==1, previousU=P.u0; else, previousU=u(k-1); end
        started = tic;

        % (4.97): kazda rzedna jest suma wazona lokalnych rzednych skoku.
        weights = ts_memberships(y(k),centers);
        s = localSteps*weights.';
        [M,freeMatrix] = dmc_matrices(s,ctrl.N,Nu);
        if nonlinearFree
            % (4.108)-(4.110) w odchyleniach od lokalnych rownowag.
            localNow = model.heq(2,:)+gain.*(previousU-equilibriumU)+pastDu.'*historyGain;
            pastWeights = ts_memberships(y(max(1,k-1)),centers);
            bias = y(k)-pastWeights*localNow.';
            for i = 1:model.count
                localFuture(:,i) = localNow(i)+localFree(:,:,i)*pastDu;
            end
            % (4.111)-(4.117): nastepna waga zalezy od poprzedniej prognozy.
            free = zeros(ctrl.N,1);
            w = weights;
            for p = 1:ctrl.N
                free(p) = localFuture(p,:)*w.'+bias;
                if p<ctrl.N, w = ts_memberships(free(p),centers); end
            end
        else
            free = y(k)+freeMatrix*pastDu;
        end

        H = 2*(M'*M+ctrl.lambda*eye(Nu));
        % Filtr przyczynowy. Regulator zna tylko biezaca wartosc zadana.
        target = alpha*target+(1-alpha)*reference(k);
        f = -2*M'*(target-free);
        bineq = [(ctrl.uMax-previousU)*ones(Nu,1); ...
                 (previousU-ctrl.uMin)*ones(Nu,1)];
        [moves,~,exitflag] = quadprog(H,f,Aineq,bineq,[],[], ...
                                      lb,ub,zeros(Nu,1),options);
        if exitflag<=0
            error('quadprog nie rozwiazal probki k=%d (exitflag=%d).',k,exitflag);
        end
        if any(Aineq*moves>bineq+1e-8) || any(moves<lb-1e-8 | moves>ub+1e-8)
            error('Plan sterowania narusza ograniczenia w probce k=%d.',k);
        end
        decisionTime(k) = toc(started);

        du = moves(1);
        u(k) = previousU+du;
        pastDu = [du;pastDu(1:end-1)];
        if k>P.delay, actualU=u(k-P.delay); else, actualU=P.u0; end
        h = rk4_step(@(state) tank_rhs(P,state,actualU,disturbance(k)), ...
                     h,P.Ts,5);
        h = max(h,0.05);
        y(k+1) = h(2);
    end
end

function dh = tank_rhs(P,h,u,d)
    f2 = P.a1*sqrt(max(h(1),0));
    f3 = P.a2*sqrt(max(h(2),0));
    cap2 = 3*P.C2*max(h(2),0.05)^2;
    dh = [(u+d-f2)/P.A1;(f2-f3)/cap2];
end
