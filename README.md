# SPM-GENETIC-TSP

## Here an example of the PMX method to form new paths

Y0 = 5713642
X0 = 4627315 

crossing-point = 3

Y1 (prime tre cifre di X e ultime 4 di Y)
x1 (prime tre cifre di Y e ultime 4 di X)

Y1 = 462|3642 -> 

La prima cifra "illegale" dopo il crossing point e' il 6, che e' ripetuto. 
Cerco il 6 nell'altro genitrore (X0) e lo sostituisco con la cifra corrispondente in X0. 
In questo caso il 6 e' in seconda posizione in X0, e lo sostituiro' con l'elemento in seconda posizione di Y0, cioe' 7

Y1 = 462|3742 

Ora anche il 4 e' duplicato. Lo sotituisco con l'ememento in prima posizione di X0 cioe' 5. 

Y1 = 462|3752

Infine stessa cosa per il 2, che corrisponde ad un 1 (terza posizione in Y0).

Y1 = 462|3751. 

Proviamo ora a calcolare anche X1: 

X1 = 571|7315 -> 571|6315 -> 571|6325 -> 571|6324
 