import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

fig, ax=plt.subplots(figsize=(8,5))
(line,)=ax.plot([],[],lw=2, color='red',label='$T(z,t)$')

ax.axhline(y=50, color='blue',linestyle='--',label='50ºC (Mort cel·lular)')
ax.axhline(y=50, color='red',linestyle='--',label='80ºC (Trombosi)')

ax.set_xlim(z[0],z[-1])
ax.set_ylim(
    np.min(T_matrix) - 2, np.max(T_matrix) + 2
)  # Marges basats en la teva matriu de temperatures
ax.set_xlabel('Posició normalitzada $\\hat{x}$')
ax.set_ylabel('Temperatura $T$ (°C)')
ax.set_title('Evolució temporal de la temperatura en el teixit')
ax.legend(loc='upper right')
ax.grid(True)

# Text opcional per mostrar el pas de temps actual a la cantonada
time_text = ax.text(
    0.05,
    0.85,
    '',
    transform=ax.transAxes,
    fontsize=12,
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7),
)