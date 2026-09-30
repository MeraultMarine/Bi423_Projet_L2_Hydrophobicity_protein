# --------------------------------------------------- #
#                  Nathan COLLIGNON                   #               
#                  Marine MERAULT                     #
#                  Kevin NGO                          #
# --------------------------------------------------- #







#------------------------Import des bibliothèque necessaires----------------------------------------------#
import tkinter as tk #Interface graphique
from tkinter import filedialog, messagebox #Pour aller ouvrir un fichier sur le PC
import numpy as np #Array pour faciliter matplotlib

import matplotlib.pyplot as plt #Graphique
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #Intégrer graphique dans l'interface tkinter



#======================================================================================#
#                                   INITIALISATION                                     #
#======================================================================================#


# Variables globales utilisées dans plusieurs fonctions
global name_tech
global graphique_dessin
global nom_fichier 

# Valeurs par défaut
graphique_dessin = False
name_tech = "Kyte-Doolittle"
nom_fichier = ""

# Dictionnaire contenant les indices d’hydrophobicité pour chaque acide aminé, pour les deux techniques
ECHELLE_HYDROPHOBE = {
    "Kyte-Doolittle": {
        "ALA":1.8,"ARG":-4.5,"ASN":-3.5,"ASP":-3.5,"CYS":2.5,
        "GLN":-3.5,"GLU":-3.5,"GLY":-0.4,"HIS":-3.2,"ILE":4.5,
        "LEU":3.8,"LYS":-3.9,"MET":1.9,"PHE":2.8,"PRO":-1.6,
        "SER":-0.8,"THR":-0.7,"TRP":-0.9,"TYR":-1.3,"VAL":4.2
    },
        "Hopp-Woods": {
        "ALA":-0.5,"ARG":3,"ASN":0.2,"ASP":3,"CYS":-1,
        "GLN":3,"GLU":0.2,"GLY":0,"HIS":-0.5,"ILE":-1.8,
        "LEU":-1.8,"LYS":3,"MET":-1.3,"PHE":-2.5,"PRO":0,
        "SER":0.3,"THR":-0.4,"TRP":-3.4,"TYR":-2.3,"VAL":-1.5
    }
}




#======================================================================================#
#                               LECTURE_PDB                                            #
#======================================================================================#


def lire_sequence_pdb(nom_f):
    '''
    Lit un fichier PDB et extrait la séquence d’acides aminés à partir des lignes SEQRES.
    Retourne une liste d’acides aminés (format 3 lettres).
    '''
    
    indice = ECHELLE_HYDROPHOBE[name_tech] #sélection de la technique demandé par l'utilisateur
    seq = []
    
    g = open(nom_f, "r")                #Ouverture du fichier en paramétre
    ligne = g.readline()  
    
    while ligne != "":
        if ligne[0:6] == "SEQRES":      #Recherche des dection SEQRES: contenant presque uniquement les acides aminées de la protéine
            long_seq = ligne[19:]       # Les AA commencent à partir du caractère 19 pour les lignes commençant par SEQRES
            i = 0
            
            while i < len(long_seq):
                aa = long_seq[i:i+3]    # Lecture par blocs de 3 lettres
                i += 4                  #saute l'acide aminé + l'espace
                if aa in indice:        #Vérifie qu'il correspond bien à un acide aminé connu
                    seq += [aa]         #Rajout à seq
        ligne = g.readline()
    g.close()
    return seq





#======================================================================================#
#                           CALCUL DU PROFIL D’HYDROPHOBICITÉ                          #
#======================================================================================#


def moyenne_hydrophobe(seq, fenetre):
    '''Prend une taille de fenêtre, un dictionnaire contenant les acides aminés et leur indice d'hydrophobicité ainsi qu'une séquence d'acides aminés pour :
    - Pour chaque acide aminé : regarder autour de celui-ci (en fonction de la taille de la fenêtre) et faire la moyenne des indices des acides aminés autour
    - Créer une liste de coordonnées x,y qu'il renverra
    - Profiter pour compter les acides aminés et rajouter un 0 à chaque fois afin de tracer à la fin la ligne de délimitation'''

    y= []                                       #listes des coordonées
    indice = ECHELLE_HYDROPHOBE[name_tech]      #sélection de la technique demandé par l'utilisateur
    
    demi_f = (fenetre-1)//2 
    position_central = demi_f                   #un coté de la fenetre (gauche)
    
    hist = False
    while position_central < len(seq)-demi_f:   # coté droit, pour chaque acide aminé. Du moment qu'il y a suffisement de place pour la fenetre (évite un out of range quand la fenetre arrive à la fin de la séquence)
        i = 0


        #########--HISTIDINE_RECHERCHE--##############
        #---> Supprimer le "clear" dans afficher_profil_hydrophobicite() si utilisation de cette partie
        #if seq[position_central] == "HIS":
        #    hist =True
        #    print("histidine à la position", position_central)
        #    zone_graphe.scatter(position_central, 0, color="red", s=10) #Verifier où est l'histidine (enlever le clear)
        ##############################################


        moyenne = 0
        while i < fenetre and position_central-demi_f+i < len(seq) : 

            if seq[position_central-demi_f+i] not in indice:            #Vérifie que l'aa existe
                i +=1
            else:
                moyenne += indice.get(seq[position_central-demi_f+i])   #va additionner les aa à droite et à gauche de notre aa principal (si fenentre = 9, on aura 4aa à droite et 4 aa à gauche)

                i +=1
        moyenne = moyenne / fenetre                                     #Calcul la moyenne des valeurs des aa


        #########--HISTIDINE_RECHERCHE--#################
        #if hist:
        #    zone_graphe.scatter(position_central, moyenne, color="red", s=10) #TESTEEEEEEE
        #    hist = False
        #################################################


        #création d'une liste comprenant toutes les coordonées des aa avec leur moyenne d'indice d'hydrophobicité
        y += [moyenne]                  #Liste des valeurs moyennes pour chaque aa
        position_central += 1           #Deplacement d'un cran en un cran à chaque fois
    
    return y                            #Retourne la listes des ordonnées pour le graphique





#======================================================================================#
#                                 GESTION DES ERREURS                                  #
#======================================================================================#


def erreur_input(taille_fenetre, sequence):
    """
    Vérifie la validité de la taille de fenêtre entrée par l’utilisateur.
    Retourne True si une erreur est détectée.
    """

    if len(sequence) == 0:
        messagebox.showerror("Erreur 01", "Aucune séquence trouvée")
        return True
    
    v = 0
    while v < len(str(taille_fenetre)):
        if str(taille_fenetre)[v] not in ["1","2","3","4","5","6","7","8","9","0"]: #convertit la taille de la fenetre en texte et vérifie que la fenêtre contient uniquement des chiffres
            messagebox.showerror("Erreur 02", "caractère non autorisé dans la taille de fenetre")
            return True
        v += 1
        
    if taille_fenetre == "":
        messagebox.showerror("Erreur 03", "Rentrez avant une taille de fenetre")
        return True
    
    elif int(taille_fenetre) == 0 or int(taille_fenetre) > len(sequence):
        messagebox.showerror("Erreur 04", "Taille de la fenetre non conforme")
        return True
    
    return False







#======================================================================================#
#                                      TUTORIEL                                         #
#======================================================================================#

def tutoriel():
    """Affiche un tutoriel d’utilisation au lancement du programme."""
    Choix = tk.messagebox.askyesno( title = "Tuto" , message = "Voulez vous un court tutoriel?" ) #Si "oui" Choix = True sinon Choix = False
    if Choix:
        messagebox.showinfo(title = "Introduction", message="Bienvenue dans ... pour vous aider à prédire les régions hydrophobes d'une protéine à partir d'un fichier PDB !")
        messagebox.showinfo(title = "Introduction", message="Pour commencer, choisissez une technique (Kyte-Doolittle ou Hopp-Woods)")
        messagebox.showinfo(title = "Introduction", message="Puis choisissez la taille de la fenêtre")
        messagebox.showinfo(title = "Introduction", message="Enfin choisissez un fichier PDB et le tour est joué !")
    
        Choix1 = tk.messagebox.askyesno( title = "interprétation" , message = "Voulez vous savoir rapidement comment interpréter le résultat ?" )
        if Choix1:
            messagebox.showinfo(title = "interprétation",message="Très bien, l'interprétation sera différente entre les deux techniques proposées.")
            messagebox.showinfo(title = "interprétation Kyte and Doolittle",message="Pour Kyte and Doolittle, les régions hydrophobes seront au dessus de la ligne, et en dessous les régions hydrophiles. Pour detecter les régions transmembranaires, une fenêtre assez large est recommandée (environ 20)")
            messagebox.showinfo(title = "interprétation Hopp and Woods",message="Pour Hopp and Woods, c'est l'inverse, les régions hydrophobes se trouvent en dessous de la ligne. Pour detecter des régions antigéniques une taille de fenetre d'environ 7 est recommandée.")





#======================================================================================#
#                                 AFFICHAGE DU GRAPHIQUE                               #
#======================================================================================#


def afficher_profil_hydrophobicite(profil, fenetre):
    '''Prend en parametre la liste des valeurs des moyenne de chaque aa
    Affiche le profil d’hydrophobicité dans la zone graphique Tkinter.'''

    global graphique_dessin         #récupere la variable global
    zone_graphe.clear()             #Efface ce qu'il y a dans la zone au cas où
    
    demi_f = (fenetre-1)//2
    positions_seq = list(range(demi_f+1, len(profil) + 1 + demi_f)) #Créer un liste de la taille de la séquence (1,2,3 etc...) en commençant par l'acide aminé à qui on a calculé sa première valeur

    valeurs = np.array(profil)

#TRACE SUR LE GRAPHIQUE
    # Courbe
    zone_graphe.plot(positions_seq, valeurs)    #Prend en parametres une liste de x et de y et pour chaque x qui correspond à un y tracer un point (et le relier)
    zone_graphe.axhline(0, linestyle="--")      #Faire une ligne de démarcation à x = 0

    # Coloration hydrophobe / hydrophile
    zone_graphe.fill_between(positions_seq, valeurs, 0, where=(valeurs >= 0), color  = "papayawhip",alpha=1) #Couleur entre graphique et valeur 0 (volume)
    zone_graphe.fill_between(positions_seq, valeurs, 0, where=(valeurs < 0),color = "blue", alpha=0.5)
    

    # Graduation X adaptée à la longueur de la séquence
    longueur = len(profil) + 2 * demi_f

    if longueur < 120: 
        step = 10
    elif longueur < 250:
        step = 20
    elif longueur < 400: 
        step = 50
    elif longueur < 1000: 
        step = 100
    elif longueur < 1500: 
        step = 200
    else: 
        step = 400
        
    zone_graphe.set_xticks(range(0, longueur + 1, step))
    zone_graphe.set_ylim(-5, 5)     #limite en ordonnée

    # Titres
    zone_graphe.set_title("Profil d'hydrophobicité : " + name_tech + " — " + nom_fichier)   #Titre
    zone_graphe.set_xlabel("Position dans la séquence (aa)")    #nom axe x
    zone_graphe.set_ylabel("Indice d'hydrophobicité")   #nom axe y

    canvas_graphique.draw()     #Trace
    graphique_dessin = True     #Variable utiliser pour signifier que le graphique à été dessiné




#======================================================================================#
#                               EXTRACTION DU NOM DE FICHIER                           #
#======================================================================================#

def nom(chemin):
    '''Extrait uniquement le nom du fichier à partir de son chemin complet.'''
    
    global nom_fichier  #récupère la variable global
    i = len(chemin)-1   #Commence au dernier indice
    nom_fichier = ""
    
    while chemin[i] != "/": #Récupere indice par indice jusqu'a tomber sur un /
        nom_fichier = chemin[i] + nom_fichier
        i -= 1


# -------- ACTION APRES APPUIE SUR BOUTON --------

#======================================================================================#
#                               CHARGEMENT DU FICHIER PDB                              #
#======================================================================================#
def charger_fichier_pdb():
    '''Ouvre un fichier PDB, calcule le profil et affiche le graphique.'''

    chemin_fichier = filedialog.askopenfilename(title="Choisir un fichier PDB",filetypes=[("Fichiers PDB", "*.pdb")])   #Demande à l'utilisateur de selectionner un fichier PDB
    taille_fenetre=fenetre_var.get() 
    
    if not chemin_fichier: 
        return  #Si aucun chemin s'arrete

    nom(chemin_fichier)
    sequence = lire_sequence_pdb(chemin_fichier)    #Lit la séquence du fichier spécifié

    if erreur_input(taille_fenetre, sequence):
        return  #Vérifie qu'aucune erreur existe avec la séquence ou la taille de la fenetre. Sinon il s'arrete
    
    taille_fenetre = int(taille_fenetre)    #Convertit la taille de fenetre en int

    label_information.config(text="Longueur : " + str(len(sequence)))   #Met à jour le bloc "longueur" de l'interface

    profil = moyenne_hydrophobe(sequence, taille_fenetre)   #calcul les y (moyenne de chaque aa)
    afficher_profil_hydrophobicite(profil, taille_fenetre)  #Dessin



#======================================================================================#
#                                   SAUVEGARDE                                          #
#======================================================================================#

def sauvegarder():
    '''Sauvegarde du graphique à l'emplacement du fichier'''
    if graphique_dessin:
        plt.savefig("profil d'hydrophobicité")
    else:
        messagebox.showerror("Erreur 05", "Aucun profil trouvé, assurez vous d'avoir généré un profil avant de sauvegarder")
        




#======================================================================================#
#                                   INTERFACE TKINTER                                   #
#======================================================================================#

# Fenêtre principale
fenetre_principale = tk.Tk()
fenetre_principale.geometry("900x650")
fenetre_principale.title("Profil d'hydrophobicité")
fenetre_principale.iconbitmap("logo.ico")
fenetre_principale.configure(bg="white")

# Configuration grille
fenetre_principale.grid_rowconfigure(0, weight=1)

fenetre_principale.grid_columnconfigure(0, weight=1)
fenetre_principale.grid_columnconfigure(1, weight=1)
fenetre_principale.grid_columnconfigure(2, weight=1)



#---------------------------------- ÉLÉMENTS DE L’INTERFACE ---------------------------#



#--------------- TITRE --------------#
titre = tk.Label(
    fenetre_principale,
    text="Profil d'hydrophobicité",
    font=("Arial", 18, "bold"),
    bg="white"
)#Lui donne son type (label text) et le texte dedans

titre.grid(row=1, column=0, columnspan=6, pady=15) #emplacement


#--------------- LABEL TAILLE FENETRE -----------------#

label_fenetre = tk.Label(
    fenetre_principale,
    text="Taille de fenêtre :",
    bg="white",
    font=("Arial", 11)
)#Lui donne son type (label text) et le texte dedans

label_fenetre.grid(row=2, column=0, padx=10, pady=10) #emplacement


#--------------- INPUT TAILLE FENETRE -----------------# 

fenetre_var = tk.StringVar()

fenetre_entry = tk.Entry(
    fenetre_principale,
    textvariable=fenetre_var,
    width=10
)#Lui donne son type (Entry) et ou il doit stocker le texte que lui donne l'utilisateur (fenetre_var)

fenetre_entry.grid(row=2, column=0, pady=10, columnspan = 2) #emplacement


#--------------- BOUTON CHARGEMENT --------------------#
bouton_chargement = tk.Button(
    fenetre_principale,
    text="Choisir fichier PDB",
    command=charger_fichier_pdb,
    bg="lightblue"
)#Lui donne son type (button), le texte dedans et ce qu'il doit faire quand il est cliqué (fonction charger_fichier_pdb)

bouton_chargement.grid(row=2, column=2, padx=10) #emplacement




#----------------- BOUTON SAVE ------------------------#
bouton_sauvergarde = tk.Button(
    fenetre_principale,
    text="Save",
    command=sauvegarder,
    bg="lightgrey"
) #Lui donne son type (button), le texte dedans et ce qu'il doit faire quand il est cliqué (fonction sauvegarder)

bouton_sauvergarde.grid(row=2, column=1, padx=10) #emplacement



# ----------------- LABEL LONGUEUR DE LA SEQUENCE ---------------#
label_information = tk.Label(
    fenetre_principale,
    text="Longueur (en aa): ",
    borderwidth=2,
    relief="solid",
    bg="lightgrey",
    width=20
)#Lui donne son type (label text) et le texte dedans

label_information.grid(row=3, column=0, padx=10,pady=10) #emplacement




# -----------------MENU CHOIX TECHNIQUES ENTRE  KD ET HOPP WOODS----------------#

def changement_k(): 
    '''Fonction appelée si K&D est choisi'''
    global name_tech
    name_tech = "Kyte-Doolittle"
    name_button_tech.set("Technique: Kyte-Doolittle")

def changement_h(): 
    '''Fonction appelée si Hopp-Woods est choisi'''
    global name_tech
    name_tech = "Hopp-Woods"
    name_button_tech.set("Technique: Hopp-Woods")


name_button_tech = tk.StringVar()
name_button_tech.set("Technique: Kyte-Doolittle")



# ---------------- Création du bouton -------------#
bouttonTech = tk.Menubutton(
    fenetre_principale,
    textvariable=name_button_tech,
    bg="lightgrey"
)

bouttonTech.grid(row=3, column=2, pady=10) #emplacement


#-------------- Création menu déroulant ---------------#
menuAffichage = tk.Menu(bouttonTech, tearoff=0)


#------------ Ajout des différentes techniques ---------------#
menuAffichage.add_command(label='Kyte-Doolittle', command=changement_k)
menuAffichage.add_command(label="Hopp-Woods", command=changement_h)
bouttonTech.configure(menu=menuAffichage)




# -------- PLACE LE GRAPH -------- #
# Création de la figure matplotlib et de sa zone de dessin (zone_graphe)
figure_graphique, zone_graphe = plt.subplots(figsize=(8, 5)) 

# Intégration du graphique matplotlib dans l’interface Tkinter
canvas_graphique = FigureCanvasTkAgg(figure_graphique, master=fenetre_principale)

# Placement du widget graphique dans la grille de la fenêtre
canvas_graphique.get_tk_widget().grid(row=4, column=0, columnspan=6, pady=20)


tutoriel()  #Propose un tutoriel quand le programme commence
fenetre_principale.mainloop()   #Boucle de la fenêtre
