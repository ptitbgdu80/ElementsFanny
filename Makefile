# Compilateur utilisé
CC=g++

# Les options de compilation
FLAGS =-std=c++11

# Le nom de l'exécutable
PROG = run

# Les fichiers source à compiler
SRC = main.cc polynomes.cpp matricesEF.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CC) $(FLAGS) $(SRC) -o $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG)
