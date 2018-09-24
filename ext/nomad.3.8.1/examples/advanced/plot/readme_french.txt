Explication programme java

Le main est situé dans la classe Prog.
Prog :
On attend de trouver le fichier init qui est créé par le c++, puis on crée la fenêtre ( le fichier init contient 0 ou 1 selon qu'on soit en simple ou bi-objectif).
Ensuite on attend de trouver le fichier in.txt qui contient les coordonnées du ou des nouveaux points.
Une fois qu'elles sont lues on renomme ce fichier en out.txt (pour que le c++ sache qu'il peut envoyer les suivants).
Les points lus sont ajoutés dans des vecteurs de points.
Une fois que les nouveaux points sont ajoutés on met à jour la fenêtre.

Fenetre :
La fenêtre est une JFrame, dont les settings de base ne devraient pas te poser de problèmes...
On a laissé la possibilité à l'utilisateur de spécifier une résolution initiale dans le fichier résolution. Sinon la fenetre est en 600 par 400.

Graph :
C'est un Jpanel à l'intérieur duquel on dessine tout.
Le constructeur est surchargé, selon qu'on soit au moment de la création de la fenetre (il faut dessiner un repère vide), en mode simple objectif(il faut ajouter un point), ou bi-objectif(on ajoute un ensemble de points).
Une fois que l'ajout a été effectué on se sert de la méthode calcul qui calcule le système de coordonnées (selon les points présents).
Puis on re dessine la fenetre (pour que le ou les nouveaux points soient présents et que le nouveau système de coordonnée soit mis en place).
Dans le cas du bi-objectif, il faut également trouver les points "extrêmes", c'est-à-dire les points dont les ordonnées ou abscisses sont les plus grande/petite, afin de connaitre la distance que notre système de coordonnées doit pouvoir couvrir.

La méthode calcul : longue et fastidieuse. 
La méthode arrondi : une méthode inutile, mais je n'ai pas trouvé comment faire autrement. Elle sert uniquement lors de l'affichage des coordonnées le long des axes. Si je n'arrondis pas les nombres avec cette méthode il ne s'affichent pas correctement (1,9997 au lieu de 2.0 par exemple).

PaintComponent :
C'est lui qui gère tout l'affichage du JPanel, et donc de la fenêtre.
Il trace les axes, leurs noms(sous forme d'images), les graduations des axes, les flèches au bout des axes, les valeurs des graduations, les points, les zones permettant l'affichage de la petite boite d'information, et la fonction en escalier.
Cette méthode est appelé lors de chaque dessin de la fenêtre.

Zone :
Permet de repérer l'arrivée de la souris dans la zone, et son départ. Lors de son arrivée elle déclenche l'affichage de la petite boite de dialogue indiquant les coordonnées du point, et sa coloration en rouge.

Voilà, c'est un peu rapide, j'espère que ça te suffira :o)

quentin
