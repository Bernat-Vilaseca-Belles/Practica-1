Com actuar amb vs code amb el github: <br>
Primer. obrir repositori. Comprovar que s'està a la carpeta correcta. Després, incertir al terminal:

    git pull

Això afegeix els nous fitxers.<br>
Després, editar el que s'hagi d'editar. <br>
Un cop s'hagi acabat d'editar els nous fitxers i s'hagin desat, s'escriu:

    git add .
    git commit -m "Add new files"
    git push

<<<<<<< HEAD
Això penja els arxius que s'ha fetgit
=======
Això penja els arxius que s'ha fet.......
>>>>>>> fdf70353f570e1dd498448cf31c4edc7a8ff6363
