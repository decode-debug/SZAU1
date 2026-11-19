function state = rk4_step(rhs,state,T,n)
    dt = T/n;
    for j = 1:n
        a = rhs(state);
        b = rhs(state+dt*a/2);
        c = rhs(state+dt*b/2);
        d = rhs(state+dt*c);
        state = state+dt*(a+2*b+2*c+d)/6;
    end
end
