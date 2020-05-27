%m-file con il quale vengono assegnati dei valori ai parametri. Questa
%operazione deve essere eseguita prima di simulare il sistema in simulink
%altrimenti i blocchi contenenti le matrici non hanno nessun valore
%definito

m1 = 2;
m2 = 2;
k = 4;
c = 5;
value = 2;
A = [0 0 1 0;0 0 0 1;-k/m1 k/m1 -c/m1 c/m1;k/m2 -k/m2 c/m2 -c/m2];
B = [0;0;1/m1;0];
C = [1 0 0 0];

[mat_S,V] = CreazioneOsservatore(A,B,C,value);