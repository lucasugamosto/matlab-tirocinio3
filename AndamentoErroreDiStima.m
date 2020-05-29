function AndamentoErroreDiStima(mat_S,x0)
    %funzione che calcola la matrice esponenziale e^(A-VC) che serve per
    %calcolare l'errore di stima durante l'intervallo interessato
    
    %presumendo che i valori della stima iniziale siano pari a 0 allora
    %l'errore di stima al tempo 0 è pari a e(0) = x(0) e inoltre sappiamo
    %che e(t) = e^((A-VC)*t)*e(0)
    
    syms s t;
    
    dim = size(mat_S);
    n = dim(1);
    
    I = eye(n);
    time = [0:1:10];
    
    mat1 = s*I;
    mat2 = mat1-mat_S;
    mat3 = inv(mat2);
    mat = ilaplace(mat3);
    
    e = mat*x0;
   
    e1 = [];
    e2 = [];
    e3 = [];
    e4 = [];
    
    for i = 0:1:10
        val = subs(e,t,i);
        
        e1 = horzcat(e1,val(1));
        e2 = horzcat(e2,val(2));
        e3 = horzcat(e3,val(3));
        e4 = horzcat(e4,val(4));
    end
    
    plot(time,e1,"b-");
    hold on
    plot(time,e2,"g-");
    hold on
    plot(time,e3,"y-");
    hold on
    plot(time,e4,"r-");
    
end