%m-file con il quale vengono assegnati dei valori ai parametri. Questa
%operazione deve essere eseguita prima di simulare il sistema in simulink
%altrimenti i blocchi contenenti le matrici non hanno nessun valore
%definito
%IN QUESTO M-FILE I VALORI NUMERICI PER I PARAMETRI (m1,m2,k,c) PERMETTONO
%DI SODDISFARE IL CASO: due autovalori immaginari puri.
fprintf("CASO: 2 autovalori immaginari puri\n");

m1 = 1;
m2 = 1;
k = 2;
c = 0;

x0 = [0;10;0;0];        %x1: posizione iniziale della massa m1
                        %x2: posizione iniziale della massa m2 (x1+L)
                        %x3: velocità iniziale della massa m1
                        %x4: velocità iniziale della massa m2

A = [0 0 1 0;0 0 0 1;-k/m1 k/m1 -c/m1 c/m1;k/m2 -k/m2 c/m2 -c/m2];
fprintf("autovalori della matrice dinamica A:\n");
eig(A)

B = [0;0;1/m1;0];
C = [1 0 0 0];
M = [0;0;0;1/m2];
value = 5;

fprintf("------------------------------------------------------------\n");
fprintf("progettazione dell'osservatore asintotico dello stato in cui\nsi impone una dinamica di errore con convergenza a zero più veloce\n");
fprintf("------------------------------------------------------------\n");

[mat_f,V] = CreazioneOsservatore(A,B,C,value);      %mat_f = A-VC

AndamentoErroreDiStima(mat_f,x0);