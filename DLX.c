#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>
#include <pthread.h>


typedef struct Thread_args_t{
    int ligneDebut;
    int ligneFin;
    int lignesZeros;
    int** relations;
    int nbLignes;
    int nbColonnes;
    int nbIncidences;
    int ** solutions;
    size_t nbSolutions;
    size_t tailleTabSolutions;
    int nbElementsParSolution;
    char* fichierEntree;
} Thread_args_t;


// structure représentant chaque élément de la matrice
typedef struct Point {
  struct Point *top;    // voisin du dessus
  struct Point *bot;    // voisin du dessous
  struct Point *left;   // voisin de gauche
  struct Point *right;  // voisin de droite
  struct Point *header; // point faisant office de tête de colonne (pas un point de la couverture)
  char type; // P = Point, R = Root (header des autres headers), H = Header
  int name;  // Nom du point si c'est un Header ou nom de la ligne si c'est un Point
  int size;  // Nombre de 1 dans la colonne (utile seulement dans les headers)
  char def;  // D = Def, N = Not Def
  int x;
  int y;
} Point;


int NB_LIGNES; // nombre de lignes dans l'espace projectif +1 pour les headers
int NB_COLONNES; // nombre de points dans l'espace projectif

Point *solutions[100];

void ajouterSolution(Thread_args_t* args, int k){
  if (args->nbSolutions == args->tailleTabSolutions){
    size_t nouvelleTaille = args->tailleTabSolutions * 2;
    int** temp = realloc(args->solutions, nouvelleTaille * sizeof(int*));
    args->solutions = temp;
    args->tailleTabSolutions = nouvelleTaille;
  }
  args->solutions[args->nbSolutions] = malloc(args->nbElementsParSolution * sizeof(int));
  for (int i = 0; i < k; i++) {
    args->solutions[args->nbSolutions][i] = solutions[i]->name;
  }
  args->nbSolutions++;
  if (args->nbSolutions%1000000 == 0){
    printf("+1 million \n");
  }
}

int** allocate_matrix(int rows, int cols) {
  int** matrix = malloc(rows * sizeof(int*));
  if (!matrix) {
      perror("Échec de l'allocation de la matrice");
      return NULL;
  }

  for (int i = 0; i < rows; i++) {
      matrix[i] = malloc(cols * sizeof(int));
      if (!matrix[i]) {
          perror("Échec de l'allocation d'une ligne");
          for (int j = 0; j < i; j++) {
              free(matrix[j]);
          }
          free(matrix);
          return NULL;
      }
  }
  return matrix;
}

int count_zero_lines(const char *filename) {
  FILE *f = fopen(filename, "r");
  if (!f) {
      perror("Impossible d'ouvrir le fichier");
      return -1;
  }

  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  int count = 0;

  while ((read = getline(&line, &len, f)) != -1) {
    bool has_zero = false;
    char *p = line;
    char *endptr;

    while (!has_zero) {
      errno = 0;
      long val = strtol(p, &endptr, 10);
      if (p == endptr) {
        // pas d'autre entier sur cette ligne
        break;
      }
      if (errno == ERANGE) {
          // valeur hors plage long, on l'ignore
      } else if (val == 0) {
        has_zero = true;
        break;
      }
      p = endptr;
    }

    if (has_zero) {
        count++;
    }
  }

  free(line);
  fclose(f);
  return count;
}

void free_matrix(int** matrix, int rows) {
  if (!matrix) return;
  for (int i = 0; i < rows; i++) {
      free(matrix[i]);
  }
  free(matrix);
}

int** lireFichier(const char* nomFichier, int* lignes, int* incidences, int* colonnes) {
  FILE* fichier = fopen(nomFichier, "r");
  if (fichier == NULL) {
      perror("Erreur lors de l'ouverture du fichier");
      return NULL;
  }
  int compteurLignes = 0;
  char ligne[1024];
  int compteurRelations = 0;
  int inWord =0;
  int nombre, max = -1;
  while (fgets(ligne, sizeof(ligne), fichier) != NULL) {
    compteurLignes++;
    if (compteurLignes == 1){
      for (int i = 0; ligne[i] != '\0'; i++) {
        if (isspace(ligne[i])) {
            inWord = 0; // On sort d'un mot
        } else if (!inWord) {
            inWord = 1; // Début d'un nouveau mot
            compteurRelations++;
        }
      }
    }
  }
  rewind(fichier);

  while (fscanf(fichier, "%d", &nombre) == 1) {
    if (nombre > max) {
        max = nombre;
    }
  }

  rewind(fichier);
  *lignes = compteurLignes;
  *incidences = compteurRelations;
  *colonnes = max+1;
  int** matrice = (int**)malloc((*lignes) * sizeof(int*));
  for (int i = 0; i < *lignes; i++) {
      matrice[i] = (int*)malloc((*incidences) * sizeof(int));
  }
  for (int i = 0; i < *lignes; i++) {
      for (int j = 0; j < *incidences; j++) {
          fscanf(fichier, "%d", &matrice[i][j]);
      }
  }
  fclose(fichier);
  return matrice;
}


void liberer_Matrice_Fichier(int** matrice, int lignes) {
  for (int i = 0; i < lignes; i++) {
      free(matrice[i]);
  }
  free(matrice);
}

void initialiser_matrice(int** matrice,int** relations, int lignes,int incidences, int colonnes) {
  for (int j = 0; j < colonnes; j++) {
    matrice[0][j] = 1;
  }

  for (int i = 1; i < lignes; i++) {
    for (int j = 0; j < colonnes; j++) {
      matrice[i][j] = 0;
    }
  }

  for (int i = 0; i < lignes-1; i++) {
    for (int j = 0; j < incidences; j++) {
      matrice[i + 1][relations[i][j]] = 1;
    }
  }
}

// création des structs Point pour chaque élément de la matrice
Point **creer_matrice(int lignes, int colonnes) {
  Point **matrice = malloc(lignes * sizeof(Point *));
  if (!matrice) {
    perror("Erreur d'allocation de la matrice");
    exit(EXIT_FAILURE);
  }
  matrice[0] = malloc((colonnes + 1) * sizeof(Point)); // ligne des headers avec 1 Header Root en plus
  if (!matrice[0]) {
    perror("Erreur d'allocation d'une ligne");
    exit(EXIT_FAILURE);
  }
  for (int i = 1; i < lignes; i++) {
    matrice[i] = malloc(colonnes * sizeof(Point));
    if (!matrice[i]) {
      perror("Erreur d'allocation d'une ligne");
      exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < lignes; i++) {
    for (int j = 0; j < colonnes; j++) {
      // on associe les voisins en prenant en compte les bords de la matrice
      // (liens circulaires)
      matrice[i][j].top =
          (i > 0) ? &matrice[i - 1][j] : &matrice[lignes - 1][j];
      matrice[i][j].bot =
          (i < lignes - 1) ? &matrice[i + 1][j] : &matrice[0][j];
      matrice[i][j].left =
          (j > 0) ? &matrice[i][j - 1] : &matrice[i][colonnes - 1];
      matrice[i][j].right =
          (j < colonnes - 1) ? &matrice[i][j + 1] : &matrice[i][0];
      matrice[i][j].header = &matrice[0][j];
      if (i == 0) {               // si c'est un header
        matrice[i][j].type = 'H'; // H = Header
        matrice[i][j].size = 0; // taille = 0 au départ puis incrémenter plus tard
        matrice[i][j].name = j;   // correspond au Point j
      } else {                    // si c'est un point
        matrice[i][j].type = 'P'; // P = Point
        matrice[i][j].name = i-1;   // relation d'incidence de la i-ème ligne
      }
      matrice[i][j].x = i;
      matrice[i][j].y = j;
      matrice[i][j].def = 'D';
    }
  }

  // initialisation du header root (il se trouve à droite du dernier header)
  matrice[0][colonnes].left = &matrice[0][colonnes - 1]; // dernier header
  matrice[0][colonnes].right = &matrice[0][0];           // premier header
  matrice[0][0].left = &matrice[0][colonnes];
  matrice[0][colonnes - 1].right = &matrice[0][colonnes];
  matrice[0][colonnes].size = 1;   // arbitraire
  matrice[0][colonnes].type = 'R'; // R = Root
  matrice[0][colonnes].x = 0;
  matrice[0][colonnes].y = colonnes;
  matrice[0][colonnes].def = 'D';

  return matrice;
}

// permet de supprimer les 0 de la matrice de points
void initialiser_liens(Point **points, int** matrice, int lignes, int colonnes) {
  for (int i = 0; i < lignes; i++) {
    for (int j = 0; j < colonnes; j++) {
      if (matrice[i][j] == 0) {
        points[i][j].left->right =
            points[i][j].right; // le voisin de droite du voisin de gauche
                                // devient le voisin de droite
        points[i][j].right->left =
            points[i][j].left; // le voisin de gauche du voisin de droite
                               // devient le voisin de gauche
        points[i][j].top->bot =
            points[i][j].bot; // le voisin du dessous du voisin du dessus
                              // devient le voisin du dessous
        points[i][j].bot->top =
            points[i][j].top; // le voisin du dessus du voisin du dessous
                              // devient le voisin du dessus
        points[i][j].size =
            0; // arbitraire mais permet d'éviter l'affichage d'un tel point
               // dans la fonction afficherVoisin()
         points[i][j].def = 'N';
      } 
      else {    // si c'est un 1
        if (i != 0) {
          points[0][j].size += 1; // taille du header de la colonne incrémentée
          points[i][j].size = 1;  // arbiraire mais permet l'affichage d'un tel
                                  // point dans la fonction afficherVoisin()
        }
      }
    }
  }
}


void insererSolutionsFichier(char* fichierSolutions, int nbThreads, Thread_args_t** args) {
  FILE *ecraserFichier = fopen(fichierSolutions, "w");
  fclose(ecraserFichier);
  FILE *fichier = fopen(fichierSolutions, "a");
    
  if (fichier == NULL) {
      perror("Erreur lors de l'ouverture du fichier !\n");
  }

  for(int i = 0 ; i < nbThreads; i++){
    for(int j = 0; j < args[i]->nbSolutions; j++){
      for(int k = 0; k < args[i]-> nbElementsParSolution; k++){
        fprintf(fichier, "%d ", args[i]->solutions[j][k]);
      }
      fprintf(fichier, "\n");
    }
  }

  fclose(fichier);

}

void liberer_matrice(Point **matrice, int lignes) {
  for (int i = 0; i < lignes; i++) {
    free(matrice[i]);
  }
  free(matrice);
}

// permet de masquer la colonne dont l'en-tête est passée en paramètres et les
// points relatifs
void masquerColonne(Point *header) {
  header->left->right = header->right;
  header->right->left = header->left;
  Point *current_bot = header->bot;
  while (current_bot != header) {
    Point *current_right = current_bot->right;
    while (current_right != current_bot) {
      current_right->top->bot = current_right->bot;
      current_right->bot->top = current_right->top;
      current_right->header->size -= 1;
      current_right = current_right->right;
    }
    current_bot = current_bot->bot;
  }
}

// permet de démasquer la colonne dont l'en-tête est passée en paramètres et les
// points relatifs
void demasquerColonne(Point *header) {
  Point *current_top = header->top;
  while (current_top != header) {
    Point *current_left = current_top->left;
    while (current_left != current_top) {
      current_left->top->bot = current_left;
      current_left->bot->top = current_left;
      current_left->header->size += 1;
      current_left = current_left->left;
    }
    current_top = current_top->top;
  }
  header->left->right = header;
  header->right->left = header;
}

void masquerLigne(int ligneI, Point **points, int colonnes){
  for (int j = 0; j < colonnes; j++) {
    if (points[ligneI][j].def == 'D') {
      points[ligneI][j].top->bot =
      points[ligneI][j].bot; // le voisin du dessous du voisin du dessus
                              // devient le voisin du dessous
      points[ligneI][j].bot->top =
      points[ligneI][j].top; // le voisin du dessus du voisin du dessous
                              // devient le voisin du dessus
      points[0][j].size -= 1;                        
    }
  }
}

void masquerLignesThread(int ligneDebut, int ligneFin, Point **points, int nbLignes, int colonnes){
  int compteurZeros = 0;
  for (int i = 1; i < nbLignes ; i++){
    if (points[i][0].def == 'D'){
      compteurZeros++;
      if (compteurZeros < ligneDebut || compteurZeros > ligneFin){
        masquerLigne(i, points, colonnes);
      }
    }
  } 
}

int DLX(Point **matrice, int k, int colonnes, Thread_args_t *args) {
  if (matrice[0][colonnes].right == &matrice[0][colonnes]) {
    ajouterSolution(args,k);
  } else {
    Point *current_colonne = matrice[0][colonnes].right;
    Point *current_choice = current_colonne;
    int max = current_choice->size;
    while (current_colonne != &matrice[0][colonnes]) {
      if (current_colonne->size < max) {
        max = current_colonne->size;
        current_choice = current_colonne;
      }
      current_colonne = current_colonne->right;
    }
    masquerColonne(current_choice);
    Point *current_point = current_choice->bot;
    while (current_point != current_choice) {
      solutions[k] = current_point;
      Point *current_point_on_line = current_point->right;
      while (current_point_on_line != current_point) {
        masquerColonne(current_point_on_line->header);
        current_point_on_line = current_point_on_line->right;
      }
      DLX(matrice, k + 1, colonnes, args);
      current_point_on_line = current_point->left;
      while (current_point_on_line != current_point) {
        demasquerColonne(current_point_on_line->header);
        current_point_on_line = current_point_on_line->left;
      }
      current_point = current_point->bot;
    }
    demasquerColonne(current_choice);
  }
}

void appel_DLX(char* fichierEntree, Thread_args_t *args)
{
  int lignes, colonnes, incidences;
  int** relations = lireFichier(fichierEntree, &lignes, &incidences, &colonnes);
  lignes +=1;
  NB_LIGNES = lignes;
  NB_COLONNES = colonnes;
  int** matrice_initiale = allocate_matrix(lignes, colonnes);
  initialiser_matrice(matrice_initiale, relations, lignes, incidences,colonnes);
  Point **matrice = creer_matrice(lignes, colonnes);
  initialiser_liens(matrice, matrice_initiale, lignes, colonnes);
  masquerLignesThread(args->ligneDebut, args->ligneFin, matrice, args->nbLignes, args->nbColonnes);
  DLX(matrice, 0, colonnes, args);
  liberer_matrice(matrice, NB_LIGNES);
  liberer_Matrice_Fichier(relations, lignes-1);
  free_matrix(matrice_initiale,lignes);
}

void *worker(void *arg) {
    Thread_args_t *args = (Thread_args_t*)arg;
    appel_DLX(args->fichierEntree, args);
    free(arg);
    return NULL;
}

int main(int argc, char** argv) {
  double time1 = (double) clock();
  time1 = time1 / CLOCKS_PER_SEC;
  char* fichierEntree = argv[1];
  char* fichierSortie = argv[2];
  int lignes, colonnes, incidences;
  
  int** relations = lireFichier(fichierEntree, &lignes, &incidences, &colonnes);
  int nbZeros = count_zero_lines(argv[1]);
  int nbThreads = atoi(argv[3]);
  int nbLignesParThread = ceil((double)nbZeros/(double)nbThreads);
  lignes +=1; // les en-têtes
  
  pthread_t threads[nbThreads];
  Thread_args_t* args[nbThreads];
  for (int i = 0; i < nbThreads; i++) {
      args[i] = malloc(sizeof(Thread_args_t));
      args[i]->ligneDebut = i*nbLignesParThread + 1;
      args[i]->ligneFin = (int)fmin((i+1)*nbLignesParThread, nbZeros);
      args[i]->lignesZeros = nbZeros;
      args[i]->nbLignes = lignes;
      args[i]->nbColonnes = colonnes;
      args[i]->nbIncidences = incidences;
      args[i]->nbElementsParSolution = args[i]->nbColonnes / args[i]->nbIncidences;
      args[i]->nbSolutions = 0;
      args[i]->tailleTabSolutions = 2;
      args[i]->solutions = malloc(args[i]->tailleTabSolutions * sizeof(int*));
      args[i]->fichierEntree = fichierEntree;
      int rc = pthread_create(&threads[i], NULL, worker, args[i]);
      if (rc) {
          fprintf(stderr, "pthread_create erreur: %d\n", rc);
          free(args[i]);
          exit(1);
      }
  }
  int nbSolutionsTotal = 0;
  for (int i = 0; i < nbThreads; i++) {
        pthread_join(threads[i], NULL);
        nbSolutionsTotal += args[i]-> nbSolutions;
  }

  insererSolutionsFichier(fichierSortie, nbThreads, args);
  sleep(2);
  double timedif = ( ((double) clock()) / CLOCKS_PER_SEC) - time1;
  printf("Le temps écoulé est de %lf secondes\n", timedif);
  printf("Nombre total de solution : %d \n", nbSolutionsTotal);
  return 0;
}
