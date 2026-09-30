# Profil d'hydrophobicité d'une protéine

Application en Python avec interface graphique (**Tkinter**) qui calcule et trace le **profil d'hydrophobicité** d'une protéine à partir d'un fichier **PDB**, afin d'aider à repérer ses régions hydrophobes (par exemple des régions transmembranaires ou antigéniques).

> Projet universitaire réalisé par **Nathan COLLIGNON**, **Marine MERAULT** et **Kevin NGO**.

## Fonctionnalités

- Lecture de la séquence d'acides aminés d'un fichier PDB (lignes SEQRES)
- Choix entre deux échelles :
  - **Kyte-Doolittle**
  - **Hopp-Woods**
- Taille de fenêtre glissante réglable par l'utilisateur
- Tracé du profil avec matplotlib intégré dans la fenêtre Tkinter :
  - ligne de référence à 0
  - zones positives et négatives colorées différemment
  - graduation de l'axe X adaptée à la longueur de la séquence
- Affichage de la longueur de la séquence
- Sauvegarde du graphique en image
- Tutoriel proposé au lancement
- Gestion des erreurs (aucune séquence trouvée, caractère non autorisé, fenêtre vide ou trop grande, sauvegarde sans profil)

## Prérequis

- Python 3.8 ou plus récent
- Bibliothèques : numpy, matplotlib (Tkinter est inclus avec Python)
- Un fichier **logo.ico** dans le même dossier que le script (utilisé comme icône de la fenêtre)


## Utilisation

1. Accepte (ou non) le petit tutoriel au lancement.
2. Choisis la **technique** dans le menu déroulant (Kyte-Doolittle par défaut).
3. Entre la **taille de la fenêtre** (nombre entier).
4. Clique sur **« Choisir fichier PDB »** et sélectionne un fichier .pdb.
5. Le profil s'affiche. Clique sur **« Save »** pour l'enregistrer.

Le graphique est enregistré sous le nom profil d'hydrophobicité.png dans le dossier de travail.

## Interpréter le résultat

| Technique | Régions hydrophobes | Taille de fenêtre conseillée |
|-----------|---------------------|------------------------------|
| **Kyte-Doolittle** | au-dessus de la ligne | ~20 (régions transmembranaires) |
| **Hopp-Woods** | en dessous de la ligne | ~7 (régions antigéniques) |

## Principe

1. Le script lit les lignes SEQRES du fichier PDB et récupère les acides aminés (code à 3 lettres).
2. Pour chaque position, il fait la **moyenne des indices** des acides aminés situés dans la fenêtre centrée sur cette position.
3. Les moyennes sont tracées en fonction de la position dans la séquence.

Les positions proches des extrémités, où la fenêtre ne tient pas entièrement, ne sont pas calculées.

## Structure du code

| Fonction | Rôle |
|----------|------|
| lire_sequence_pdb() | Extrait la séquence depuis les lignes SEQRES |
| moyenne_hydrophobe() | Calcule la moyenne glissante des indices |
| erreur_input() | Vérifie la séquence et la taille de fenêtre |
| tutoriel() | Affiche le tutoriel de démarrage |
| afficher_profil_hydrophobicite() | Trace le graphique dans l'interface |
| nom() | Extrait le nom du fichier depuis son chemin |
| charger_fichier_pdb() | Enchaîne lecture, calcul et affichage |
| sauvegarder() | Enregistre le graphique |

## Références

- Kyte J. & Doolittle R.F. (1982), *A simple method for displaying the hydropathic character of a protein*, J. Mol. Biol.
- Hopp T.P. & Woods K.R. (1981), *Prediction of protein antigenic determinants from amino acid sequences*, PNAS.
- [Protein Data Bank (PDB)](https://www.rcsb.org/)
