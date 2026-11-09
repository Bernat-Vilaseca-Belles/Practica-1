#Càlculs per realitzar el mètode d'Euler Explícit.

Tenim:
$$\partial_{\tilde t} \tilde T= \partial^2_{\tilde z} \tilde T + 1$$
Tenim una discretització de $\tilde z$ de $N=101$ punts. Per tant:
$$\tilde z = \frac z {L_0} \Rightarrow \Delta \tilde z = \frac{\Delta z}{L_0} \xrightarrow{\Delta z = \frac{h-0}{N-1}} \Delta \tilde z= \frac h h \frac{1}{N-1}= \frac 1{100}$$
Tenim discretitzacions en el temps de:
    > $\delta \tilde t= 0.51 (\delta \tilde z)^2, \quad \delta \tilde t= 0.49 (\delta \tilde z)^2, \quad \delta \tilde t= 0.51 (\delta \tilde z)^2$
Enomenem $\phi_j^n$ a la temperatura normalitzada a temps $n$ i posició en $z$ $j$.
Volem trobar una expressió que ens doni el valor de $T$ a temps $n+1$.
$$\partial_{\tilde t} \tilde T(\tilde z_j, \tilde t_n)\approx \dot \phi_j^n = \frac{\phi_j^{n+1}- \phi_j^n}{\Delta \tilde t}$$

$$\partial_{\tilde z}^2 \tilde T(\tilde z_j, \tilde t_n)\approx ((\phi_j^n)')'$$
$$(\phi_j^n)' = \frac{\phi_{j+1}^n - \phi_j^n}{\Delta \tilde z} \quad \text{(Derivada per la dreta)}$$
$$((\phi_j^n)')'= \frac{(\phi_{j}^n)' - (\phi_{j-1}^n)'}{\Delta \tilde z} \quad \text{(Derivada per l'esquerra)}$$
$$((\phi_j^n)')'= \frac{\phi_{j+1}^n - \phi_j^n- \phi_j^n+ \phi_{j-1}^n}{(\Delta \tilde z)^2}=\frac{\phi_{j+1}^n - 2\phi_j^n+ \phi_{j-1}^n}{(\Delta \tilde z)^2}$$
Ho apliquem a l'equació:
$$\frac{\phi_j^{n+1}- \phi_j^n}{\Delta \tilde t}= \frac{\phi_{j+1}^n - 2\phi_j^n+ \phi_{j-1}^n + (\Delta \tilde z)^2}{(\Delta \tilde z)^2}$$
$$\phi_j^{n+1}=\phi_j^n + \gamma(\phi_{j+1}^n - 2\phi_j^n+ \phi_{j-1}^n + (\Delta \tilde z)^2),$$
on $\gamma= \frac{\Delta \tilde t}{(\Delta \tilde z)^2}$. És estable sempre que $\gamma \le 1/2$. Per tant, tindrem un sistema inestable i caòtic amb la primera discretització. A més, aquests càlculs són per $n=1,...N-1$ i $j=2,..., N-1$. Ho hem calculat així perquè sabem sempre que la temperatura no normalitzada als punts extrems és $36.5^\circ \text C$.