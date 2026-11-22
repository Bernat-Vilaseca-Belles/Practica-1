# Practica-1
## Declaració d'ús de la Intel·ligència Artificial (IA):
Nosaltres declarem que només hem utilitzat la intel·ligència artificial com a font d'informació que després ha sigut contrastada amb Stack Overflow i les biblioteques d'informació de Python; com a corrector del tipus ortogràfic i com a consulta de funcions per a la redacció de l'informe a llenguatge LaTeX.

## Com s'empra aquest programa
Nosaltres recomanem que es guardi en una carpeta en format local i s'obri amb el programari Visual Studio Code. Majoritàriament no requereix d'un ordre en específic d'execució però sí que hi ha uns pocs programes que s'han d'executar en aquest ordre:
1. Euler_explícit.py
2. Euler_implícit.py
3. Crank-Nicolson.py
4. Solució_problema.py
5. animació1.py

## Autoria del treball
A continuació es relata qui és l'autor o els autors de cada part.
- Euler Explícit: Bernat Vilaseca.
- Gauss-Seidel: Javier Álvarez, Joel Serrano i Bernat Vilaseca.
- Euler Implícit: Javier Álvarez, Joel Serrano i Bernat Vilaseca.
- Crank-Nicolson: Javier Álvarez i Joel Serrano.
- Solució del problema: Bernat Vilaseca (emprant el procediment dut a terme pels altres integrants del grup en el mètode de Crank-Nicolson).
- Animació: Joel Serrano i Bernat Vilaseca.
- Càlculs: Javier Álvarez i Joel Serrano (amb contribució de Bernat Vilaseca).

## Què conté aquest repositori?
Els continguts es llisten a continuació:
- CN.png: El gràfic resultant d'utilitzar el mètode de Crank-Nicolson.
- Comparació EE_CN.png: El gràfic resultant de la comparació entre el mètode d'Euler Explícit (amb cada discretització) amb el Crank-Nicolson.
- Comparació EI_CN.png: El gràfic resultant de la comparació entre el mètode d'Euler Implícit (amb cada discretització) amb el Crank-Nicolson.
- Crank-Nicolson.py: El codi que implementa el mètode de Crank-Nicolson.
- EE.png: El gràfic resultant d'implementar el mètode d'Euler Explícit.
- EI.png: El gràfic resultant d'implementar el mètode d'Euler Implícit.
- Euler_explícit.py: El codi que implementa el mètode d'Euler Explícit.
- Euler_implícit.py: El codi que implementa el mètode d'Euler Implícit.
- README.md: Aquest mateix document que vostè està llegint ara mateix.
- Resultat_solució.png: El gràfic de la temperatura quan s'atura la cirurgia (quan s'acaba la simulació), implementat amb el mètode de Crank-Nicolson.
- Solució_problema.py: El codi que simula tota la cirurgia amb el mètode de Crank-Nicolson.
- Temperatura_ee.csv: Llistat de la temperatura a cada punt a $\tilde t = 0.025$ per a cada desratització utilitzant el mètode d'Euler Explícit.
- Temperatura_ei.csv: Llistat de la temperatura a cada punt a $\tilde t = 0.025$ per a cada desratització utilitzant el mètode d'Euler Implícit.
- Videoanimacio.gif: El video de com evoluciona temporalment la temperatura, representat en funció del color, per una regió determiada i amb el mètode de Crank-Nicolson.
- animacio1.py: Codi que representa de manera visual el que passa en el document Solució_problema.py.
- constants.py: Document amb totes les constants que fem servir nosaltres en els altres codis.
- gauss_seidel.py: Mètode iteratiu que s'ha emprat per trobar la solució d'un sistema d'equacions per poder resoldre els mètodes d'Euler Implícit I Crank-Nicolson.
- solució_analítica.py: Document que conté una funció que retorna la temperatura no normalitzada en Kelvin de tots els punts en un temps determinat.
