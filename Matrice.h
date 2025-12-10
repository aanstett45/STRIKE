#ifndef MATRICE_H
#define MATRICE_H

#include <map>
#include <vector>
#include <string>
#include <omp.h>

#include "Papier.h"
#include <roaring/roaring.hh>

using RoaringBitmap = Roaring;
using namespace std;

class Matrice {
private:
    map<int, Papier*> matrice;
    int nb_papier;

    // Méthodes utilitaires internes
    int dot_with_restart(const vector<int>& rowIndices, const vector<int>& colIndices,
                         Matrice* csrMatrix, Matrice* cscMatrix, int rowIndex, int colIndex, int& colPtrRef);
public:
    // Constructeur par défaut
    Matrice();
    Matrice(const Matrice &m) = default;
    ~Matrice();

    // Méthodes publiques (implémentées dans Matrice.cpp)
    void AddPapier(int id, const string &titre);
    Papier *GetPapier(int id);
    void AddLigne(int val, int idArticle, int idCite);
    void AddColumns(int idArticle, int idCitation);
    void AddToLigne(int val, int idArticle, int idCite);
    bool IsInMatrice(int id);
    int GetNbPapier() const;

    int GetNnz() const;

    // Multiplications implémentées
    Matrice* Multiply(Matrice* m);
    Matrice* MultiplyNaive(Matrice* m);
    Matrice* MultiplyWithRestart(Matrice* m);
    Matrice* MultiplyParallelWithRestart(Matrice* m);

    // Addition et copie
    void Add(Matrice *m);
    void CopyAllPaper(Matrice *m);

    // Utilitaire (déclaré mais implémentation éventuelle)
    void SetDiagonalToZero();
};

#endif // MATRICE_H
#endif // MATRICE_H