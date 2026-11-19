%% Niezalezny test predykcji NPL przez splot przyrostow odpowiedzi skokowej z wejsciem
clear; clc;
load(fullfile('wyniki_czesc4','dane_czesc4.mat'),'results','ctrl','Nu','Tf');
load(fullfile('wyniki_czesc3','dane_czesc3.mat'),'P','modelTuned','localSteps','traces');
D=ctrl.D; N=ctrl.N; impulse=diff([zeros(1,modelTuned.count);localSteps]);
worst=0; opt=optimoptions('lsqlin','Display','off');
for j=1:numel(results)
    y=results(j).y; u=results(j).u; ref=traces(j).ref;
    pole=exp(-P.Ts/Tf); target=filter(1-pole,[1,-pole],ref,pole*ref(1));
    padded=[repmat(P.u0,D,1);u]; eq=modelTuned.qeq-P.d0;
    for k=1:numel(u)
        prev=padded(D+k-1);
        % p=0 to chwila pomiaru, p>0 zaklada stale przyszle sterowanie.
        indices=min(D+k+(0:N)'-(1:D),D+k-1);
        local=padded(indices)*impulse+modelTuned.heq(2,:)-eq.*localSteps(end,:);
        w=ts_memberships(y(max(1,k-1)),modelTuned.centers);
        bias=y(k)-local(1,:)*w.';
        z=y(k); free=zeros(N,1);
        for p=1:N
            w=ts_memberships(z,modelTuned.centers);
            z=local(p+1,:)*w.'+bias; free(p)=z;
        end
        s=localSteps*ts_memberships(y(k),modelTuned.centers).';
        M=toeplitz(s(1:N),[s(1),zeros(1,Nu-1)]);
        G=[M;sqrt(ctrl.lambda)*eye(Nu)]; b=[target(k)-free;zeros(Nu,1)];
        plan=G\b; L=tril(ones(Nu)); bounds=[repmat(ctrl.uMax-prev,Nu,1);repmat(prev-ctrl.uMin,Nu,1)];
        if any(abs(plan)>ctrl.duMax) || any([L;-L]*plan>bounds)
            lim=ctrl.duMax*ones(Nu,1);
            [plan,~,~,flag]=lsqlin(G,b,[L;-L],bounds,[],[],-lim,lim,[],opt); assert(flag>0);
        end
        worst=max(worst,abs(plan(1)-(u(k)-prev)));
    end
    assert(all(u>=ctrl.uMin-1e-8 & u<=ctrl.uMax+1e-8) && max(abs(diff([P.u0;u])))<=ctrl.duMax+1e-8);
end
fprintf('Najwieksza roznica ruchu NPL wzgledem niezaleznego splotu i LS: %.3g\n',worst);
assert(worst<1e-6);
