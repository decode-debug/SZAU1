%% Niezalezna weryfikacja implementacji czesci 4 (FGPC-NPL, wzory 4.118-4.125)
clear; clc;
load(fullfile('wyniki_czesc4','dane_czesc4.mat'), ...
    'nplTraces','nplValTraces','N_npl','Nu_npl','lambdaNPL','ctrl');
load(fullfile('wyniki_czesc3','dane_czesc3.mat'), ...
    'P','modelTuned','traces','tControl','refStart','distStart','validationCases');

saved = [nplTraces;nplValTraces];
refs = {traces.ref};
disturbances = {traces.dist};
for j = 1:size(validationCases,1)
    refs{end+1} = P.h0(2)+validationCases(j,1)*(tControl>=refStart); %#ok<SAGROW>
    disturbances{end+1} = P.d0+validationCases(j,2)*(tControl>=distStart); %#ok<SAGROW>
end

worstMove = 0;
worstConstraint = 0;
worstPlant = 0;
opt = optimoptions('lsqlin','Display','off');
L = tril(ones(Nu_npl));
onesU = ones(Nu_npl,1);
lim = ctrl.duMax*onesU;

for scenario = 1:size(saved,1)
    ref = refs{scenario};
    dist = disturbances{scenario};
    for v = 1:numel(lambdaNPL)
        y = saved(scenario,v).y;
        u = saved(scenario,v).u;
        du = diff([P.u0;u]);
        worstConstraint = max([worstConstraint;u-ctrl.uMax;ctrl.uMin-u;abs(du)-ctrl.duMax]);

        % Niezalezna rekonstrukcja iteracyjnej predykcji swobodnej (4.118)-(4.125) i QP
        for k = 1:numel(u)
            if k==1, prev = P.u0; else, prev = u(k-1); end
            % Macierz dynamiczna M_k z linearyzacji w y(k)
            w = ts_memberships(y(k),modelTuned.centers);
            b = modelTuned.b*w.'; c = modelTuned.c*w.';
            s = zeros(N_npl,1);
            for p = 1:N_npl
                y1 = 0; y2 = 0;
                if p>1, y1 = s(p-1); end
                if p>2, y2 = s(p-2); end
                s(p) = b(1)*y1 + b(2)*y2 + c(1)*(p>P.delay) + c(2)*(p>P.delay+1);
            end
            M = toeplitz(s,[s(1),zeros(1,Nu_npl-1)]);

            % (4.118)-(4.119): biezace wyjscie modelu rozmytego y^M_k i d_k
            ykm1 = y(max(1,k-1)); ykm2 = y(max(1,k-2));
            w0 = ts_memberships(ykm1,modelTuned.centers);
            yMk = (modelTuned.b*w0.')'*[ykm1;ykm2] + ...
                  (modelTuned.c*w0.')'*[get_u(u,k-1-P.delay,k-1,prev,P.u0); ...
                                        get_u(u,k-2-P.delay,k-1,prev,P.u0)] + ...
                  modelTuned.a0*w0.';
            dk = y(k) - yMk;

            % (4.120)-(4.125): rekurencyjna nieliniowa trajektoria swobodna
            free = zeros(N_npl,1); y1 = y(k); y2 = ykm1;
            for p = 1:N_npl
                wp = ts_memberships(y1,modelTuned.centers);
                yNext = (modelTuned.b*wp.')'*[y1;y2] + ...
                        (modelTuned.c*wp.')'*[get_u(u,k+p-1-P.delay,k-1,prev,P.u0); ...
                                              get_u(u,k+p-2-P.delay,k-1,prev,P.u0)] + ...
                        modelTuned.a0*wp.' + dk;
                free(p) = yNext; y2 = y1; y1 = yNext;
            end

            bounds = [(ctrl.uMax-prev)*onesU;(prev-ctrl.uMin)*onesU];
            if lambdaNPL(v)==0
                [plan,~,~,flag] = lsqlin(M,ref(k)-free,[L;-L],bounds,[],[],-lim,lim,[],opt);
                assert(flag>0);
            else
                G = [M;sqrt(lambdaNPL(v))*eye(Nu_npl)];
                rhs = [ref(k)-free;zeros(Nu_npl,1)];
                plan = G\rhs;
                if any(abs(plan)>lim) || any([L;-L]*plan>bounds)
                    [plan,~,~,flag] = lsqlin(G,rhs,[L;-L],bounds,[],[],-lim,lim,[],opt);
                    assert(flag>0);
                end
            end
            worstMove = max(worstMove,abs(plan(1)-du(k)));
        end

        % Niezalezny bilans objetosci dla obiektu
        vol = [P.A1*P.h0(1);P.C2*P.h0(2)^3];
        chk = zeros(numel(u)+1,1); chk(1) = P.h0(2);
        for k = 1:numel(u)
            if k>P.delay, actU = u(k-P.delay); else, actU = P.u0; end
            dt = P.Ts/5;
            for sub = 1:5
                a = vol_rhs(P,vol,actU,dist(k));
                b = vol_rhs(P,vol+dt*a/2,actU,dist(k));
                c = vol_rhs(P,vol+dt*b/2,actU,dist(k));
                d = vol_rhs(P,vol+dt*c,actU,dist(k));
                vol = vol+dt*(a+2*b+2*c+d)/6;
            end
            chk(k+1) = (vol(2)/P.C2)^(1/3);
        end
        worstPlant = max(worstPlant,max(abs(chk-y)));
    end
end

fprintf('Maks. naruszenie ograniczen (FGPC-NPL): %.3g\n',worstConstraint);
fprintf('Maks. blad ruchu QP wzgledem niezaleznego lsqlin: %.3g\n',worstMove);
fprintf('Maks. roznica symulacji wysokosci i objetosci: %.3g cm\n',worstPlant);
assert(worstConstraint<1e-8 && worstMove<1e-6 && worstPlant<2e-3);
fprintf('AUDYT CZESCI 4 (FGPC-NPL): WSZYSTKIE TESTY ZALICZONE.\n');

function val = get_u(u,idx,lastKnown,prev,u0)
    if idx<1, val = u0; elseif idx<=lastKnown, val = u(idx); else, val = prev; end
end

function dv = vol_rhs(P,v,u,d)
    h1 = max(v(1)/P.A1,0);
    h2 = max(v(2)/P.C2,0)^(1/3);
    dv = [u+d-P.a1*sqrt(h1); P.a1*sqrt(h1)-P.a2*sqrt(h2)];
end
