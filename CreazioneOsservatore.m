function [mat_f,V] = CreazioneOsservatore(A,B,C,value)
    %inizializzazione dei valori da utilizzare
    syms t
    dim_A = size(A);
    n = dim_A(1);
    I = eye(n);
    
    %controllo osservabilità del sistema primale S
    Q = [];
    for i = 0:n-1
        Q = vertcat(Q,C*(A^i));
    end
    
    if rank(Q) == n
        %fprintf("il sistema primale S è osservabile\n");
        oss = 1;
    else
        fprintf("il sistema primale S non è osservabile\n");
        oss = 0;
    end
    
    %controllo determinabilità del sistema primale S
    if oss == 1
        %fprintf("il sistema primale S è determinabile\n");
    else
        mat_A = A^4;
        mat_B = vertcat(Q,mat_A);
        
        if rank(mat_B) == rank(Q)
            %fprintf("il sistema primale S è determinabile\n");
            det = 1;
        else
            fprintf("il sistema primale S non è determinabile\n");
            det = 0;
        end
    end
    
    %controllo rilevabilità del sistema primale S
    if oss == 0 && det == 0
        fprintf("non è possibile progettare l'osservatore per S\n");
    else
        autovalori_A = eig(A);
        for i = 1:n
            if real(autovalori_A(i)) >= 0
                mat_1 = autovalori_A(i)*I;
                mat_2 = A-mat_1;
                mat_3 = vertcat(mat_2,C);
                
                if rank(mat_3) ~= n
                    fprintf("non è possibile progettare l'osservatore per S\n");
                    return
                end
            end
        end
        
        fprintf("è possibile progettare l'osservatore per il sistema S\n");
    end
                
    %progettazione dell'osservatore per il sistema primale S
    
    Ad = A';
    Bd = C';
    Cd = B';
    Pd = Q';
    
    if oss == 1
        %caso in cui il sistema primale S è osservabile
        
        %fprintf("il sistema duale Sd è raggiungibile\n");
        
        %calcolo del vettore tau
        vet = [];
        for i = 1:n
            if i == n
                vet = horzcat(vet,1);
            else
                vet = horzcat(vet,0);
            end
        end            
        tau = vet*(inv(Pd));
        
        %calcolo del polinomio caratteristico desiderato
        newAutovalori = [];
        for i = 1:n
            val = real(autovalori_A(i));
            val = round(val);
            val = val-value;
            newAutovalori = vertcat(newAutovalori,val);
        end
        
        Pdes = 1;
        for i = 1:n
            val = t-newAutovalori(i);
            Pdes = Pdes*val;
        end
        Pdes = expand(Pdes);
        vet = sym2poly(Pdes);
        
        %sostituzione della variabile t con la matrice Ad per il calcolo di
        %Fd
        mat = zeros(n);
        j = 0;
        for i = n+1:-1:1
            val = vet(i)*(Ad^j);
            mat = mat+val;
            j = j+1;
        end    
        
        %calcolo di Fd e della rispettiva V
        Fd = -tau*mat;
        V = -(Fd');

        mat_f = A-(V*C);
    else
        %caso in cui il sistema primale S non è osservabile
        
        %fprintf("il sistema duale Sd non è raggiungibile\n");
        
        %occorre preliminarmente calcolare le matrici T e T^-1 per la forma
        %di kalman sul sistema duale Sd
        
        %calcolo di una base per Xr
        Xr = [];
        dim_Xr = 0;
        nr = rank(Pd);
        j = 1;
        
        while dim_Xr < nr
            if Pd(:,j) == 0
                j = j+1;
            else
                Xr = horzcat(Xr,Pd(:,j));
                dim_Xr = dim_Xr+1;
                j = j+1;
            end
        end
        
        %calcolo del completamento di Xr
        j = 1;
        
        while dim_Xr < n
            p = horzcat(Xr,I(:,j));
            if rank(p) == n
                dim_Xr = dim_Xr+1;
                j = j+1;
            else
                p = Xr;
                j = j+1;
            end
        end
        T_inv = p;
        T = inv(T_inv);
        
        %calcolo delle nuove matrici nella forma di kalman
        newAd = T*Ad*T_inv;
        newBd = T*Bd;
        
        %calcolo delle matrici della parte raggiungibile del sistema Sd
        Adr = [];
        Bdr = [];
        
        for i = 1:nr
            Adr = horzcat(Adr,newAd(1:nr,i));
            Bdr = vertcat(Bdr,newBd(i,:));
        end
        
        %calcolo degli autovalori appartententi al sottosistema
        %raggiungibile di Sd e dei nuovi autovalori desiderati
        autovalori_Adr = eig(Adr);
        newAutovalori = [];
        
        for i = 1:nr
            val = real(autovalori_Adr(i));
            val = round(val);
            val = val-value;
            newAutovalori = vertcat(newAutovalori,val);
        end
        
        %calcolo della matrice di raggiungibilità per il sottosistema
        %raggiungibile di Sd
        Pdr = [];
        
        for i = 0:nr-1
            Pdr = horzcat(Pdr,(Adr^i)*Bdr);
        end
        
        %calcolo della matrice Fd per il sottosistema raggiungibile di Sd
        vet = [];
        
        for i = 1:nr
            if i == nr
                vet = horzcat(vet,1);
            else
                vet = horzcat(vet,0);
            end
        end            
        tau = vet*(inv(Pdr));
        
        Pdes = 1;
        for i = 1:nr
            val = t-newAutovalori(i);
            Pdes = Pdes*val;
        end
        Pdes = expand(Pdes);
        vet = sym2poly(Pdes);
        
        mat = zeros(nr);
        j = 0;
        for i = nr+1:-1:1
            val = vet(i)*(Adr^j);
            mat = mat+val;
            j = j+1;
        end
        
        %calcolo di Fd e di V
        Fdr = -tau*mat;
        vet = zeros(n-nr,1);
        
        while nr < n
            Fdr = horzcat(Fdr,vet);
            nr = nr+1;
        end
        Fd = Fdr*T;
        
        V = -(Fd');
        
        mat_f = A-(V*C);
    end
    
    fprintf("i nuovi autovalori della matrice A-VC sono:\n");
    eig(mat_f)
    %controllo esattezza del risultato trovato
    %se il polinomio caratteristico di A-VC è uguale a quello desiderato
    %allora il risultato è giusto
    autovalori_mat_f = eig(mat_f);
    P_mat_S = 1;
    for i = 1:n
        val = t-autovalori_mat_f(i);
        P_mat_S = P_mat_S*val;
    end
    P_mat_S = expand(Pdes);
    
    vet1 = sym2poly(Pdes);
    vet2 = sym2poly(P_mat_S);
    
    if vet1 == vet2
        fprintf("Pa-vc == Pdes -> la matrice A-VC è corretta\n");
    else
        fprintf("Errore, i polinomi sono diversi\n");
    end
end