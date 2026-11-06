Anomenem $h$ a la longitud del cilíndre que estem estudiant, expressada en metres. $h= 0.02 \, \text{m}$<br>
Nosaltres ho hem calculat amb dependència amb $z$, altura.
$$\partial_t T = \frac{k}{c_v \rho} \nabla^2 T + \frac{1}{c_v \rho} P_{ext}$$

**Càlcul de $P_{ext}$**
$$V=V_0 \cos{(\omega t)}$$
$$V'= V_0/\sqrt 2 \text{ (potencial efectiu)}$$
Potència dissipada per una resistència en corrent altern:
$$P=V'^2/R,$$
on $R$ és la resistència i $P$ és la potència.
$$R=\frac h {S \sigma},$$
on $\sigma$ és la conductivitat (uniforme en el nostre medi) i S és la secció.
$$P=\frac{V_0^2 S \sigma}{2h}$$
Podem deduïr fàcilment que $P_{ext}= \frac P {\text{Volum}}=\frac{V_0^2 S \sigma}{2 h^2 S}= \frac{V_0^2 \sigma}{2 h^2}$


$$\nabla^2 T= \frac{1}{r} \frac \partial {\partial r} \left(r \frac{\partial T}{\parital r} \right)+ \frac 1 {r^2} \frac{\partial^2 T}{\partial \phi^2} + \frac{\partial^2 T}{\partial z^2}$$
Sabem que no depèn del radi ni de l'angle. així que:
$$\nabla^2 T= \partial_z^2 T$$

$$\partial_t T = \frac{k}{c_v \rho} \partial_z^2 T + \frac{V_0^2\sigma}{2\rho C_v h^2}$$

**Normalització:**
*Cal definir què és normalització.*
Fem les següents anomenacions:
    - $T' = \frac T {T_0}$
    - $z' = \frac z {L_0}$
    - $t' = \frac t {t_0}$

Denotem  que $\alpha= \frac{k}{c_v \rho}$. Operant, arribem:
$$\frac{\partial T'}{\partial t'} = \textcolor{blue}{\alpha \frac{t_0}{L_0^2}} \frac{\partial^2 T'}{\partial z'^2}+ \textcolor{green}{\frac{t_0 \alpha}{t_0 k}\frac{V'}{2h^2}}$$

Aquí podem fer diverses suposicions per tal de normalitzar-ho. S'han de fer 3, i es pot escollir què suposar. Serveix per determinar el valor dels paràmetres $t_0, T_0, L_0$.<br>
Nosaltres escollim que $L_0 =h$ i que la totalitat dels termes en blau i en verd sigui 1.
Així ens queda:
$$L_0 =h; \quad T_0 = \frac k {\sigma V_0^2}; \quad t_0=\frac{h^2 c_v \rho}{k}$$
$$\boxed{\partial_{t'} T= \partial^2_{z'} T + 1}$$
Aquí veiem que el terme $q(z,t)$ de l'anex del document enunciat és constant en el temps i en l'espai.

*(Per fer aquest càlcul no hem tingut en compte potencials retardats, degut al fet que estem en quasi-estàtica. Sinó, $q(z,t)$ sí que ens dependria de l'espai i del temps)*