//
// Created by alex on 22/03/24.
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <omp.h>
#include <queue>

#include <algorithm>
#include <chrono>

#include "Matrice.h"
#include "Papier.h"

#include <roaring/roaring.hh>

#include <random>
#include <cmath>

using RoaringBitmap = Roaring;
using namespace std;

Matrice::Matrice(){
    this->nb_papier = 0;
}

int Matrice::GetNnz()const {
    int mnnz = 0;
    for( auto elem : matrice){
        mnnz += elem.second->GetNnz();
    }
    return mnnz;
}


void Matrice::AddPapier(int id, const string &titre) {
    if(IsInMatrice(id)){
        this->GetPapier(id)->SetTitre(titre);
    }else{
        this->matrice[id] = new Papier(titre);
        this->nb_papier += 1;
    }
}

Papier *Matrice::GetPapier(int id) {
    return matrice[id];
}

void Matrice::AddLigne(int val, int idArticle, int idCite) {
    this->matrice[idArticle]->AddLigne(val,idCite);
}

void Matrice::AddColumns(int idArticle,int idCitation){
    if(IsInMatrice(idCitation)){
        this->GetPapier(idCitation)->AddColumn(idArticle);
    }else{
        this->AddPapier(idCitation,"Titre indisponible");
        this->GetPapier(idCitation)->AddColumn(idArticle);
    }
}



// Gustavson naïf : utilise accum + touched (pas d'optimisation bitmap/roaring)
Matrice* Matrice::MultiplyNaive(Matrice* m) {
    auto mtmp = new Matrice();
    mtmp->CopyAllPaper(this);

    const int N = this->GetNbPapier();

    // accumulateur réutilisé pour chaque ligne
    std::vector<int> accum(N, 0);
    std::vector<int> touched;
    touched.reserve(1024);

    for (int i = 0; i < N; ++i) {
        if (i % 10000 == 0) {
            std::cout << "MultiplyNaive progress: " << (i * 100) / std::max(1, N) << " %\n";
        }

        touched.clear();

        // récupère la liste des k où A[i,k] != 0 (sans Roaring)
        const std::vector<int> rowKs = this->GetPapier(i)->GetAllLignes();
        if (rowKs.empty()) continue;

        // pour chaque k non nul dans la ligne i
        for (int k : rowKs) {
            int a_val = this->GetPapier(i)->GetValeurById(k);
            if (a_val == 0) continue;

            // récupère la liste des j où B[k,j] != 0 (sans Roaring)
            const std::vector<int> colsFromK = m->GetPapier(k)->GetAllColumns();
            for (int j : colsFromK) {
                int b_val = m->GetPapier(k)->GetValeurById(j);
                if (b_val == 0) continue;

                int old = accum[j];
                int newv = old + a_val * b_val;
                if (old == 0 && newv != 0) touched.push_back(j);
                accum[j] = newv;
            }
        }

        // écrire les résultats pour la ligne i
        for (int col_j : touched) {
            int v = accum[col_j];
            if (v != 0) {
                mtmp->AddToLigne(v, i, col_j);
            }
            accum[col_j] = 0;
        }
    }

    return mtmp;
}


Matrice* Matrice::Multiply(Matrice* m) {
    // mtmp = résultat C
    auto mtmp = new Matrice();
    mtmp->CopyAllPaper(this); // j'assume que ça initialise la structure résultat

    const int N = this->GetNbPapier(); // dimension carrée supposée

    // Accumulateur temporaire pour une ligne i de la sortie
    // IMPORTANT : on l'alloue une fois et on le réutilise pour chaque i
    std::vector<int> accum(N, 0);
    std::vector<int> touched;
    touched.reserve(1024); // petit réservoir pour limiter les reallocs

    for (int i = 0; i < N; i++) {
        if (i % 10000 == 0) {
            std::cout << "Multiply progress: " << (i * 100) / N << " %\n";
        }

        touched.clear();

        // 1. Récupère les colonnes k où A[i,k] != 0
        const Roaring& rowNZ = this->GetPapier(i)->getLignesBitmap();
        // rowNZ = indices k actifs dans la ligne i de A

        // 2. Pour chaque k non nul dans la ligne i
        for (uint32_t k32 : rowNZ) {
            int k = static_cast<int>(k32);

            // valeur a = A[i,k]
            int a_val = this->GetPapier(i)->GetValeurById(k);
            if (a_val == 0) continue; // just in case

            // 3. Propage via la "ligne k" de B:
            // getColumnsBitmap() doit donner l'ensemble des j tels que B[k,j] != 0
            const Roaring& colsFromK = m->GetPapier(k)->getColumnsBitmap();

            for (uint32_t j32 : colsFromK) {
                int j = static_cast<int>(j32);

                // valeur b = B[k,j]
                int b_val = m->GetPapier(k)->GetValeurById(j);
                if (b_val == 0) continue; // sécurité

                int old_val = accum[j];
                int new_val = old_val + a_val * b_val;
                if (old_val == 0 && new_val != 0) {
                    touched.push_back(j); // première fois qu'on touche cette colonne j pour la ligne i
                }
                accum[j] = new_val;
            }
        }

        // 4. On émet la ligne i de C à partir de accum
        for (int col_j : touched) {
            int v = accum[col_j];
            if (v != 0) {
                // C[i,j] = v
                mtmp->AddToLigne(v, i, col_j);
            }
            // reset pour réutiliser accum proprement à la ligne suivante
            accum[col_j] = 0;
        }
    }

    return mtmp;
}



int Matrice::dot_with_restart(const std::vector<int>& rowIndices,
                              const std::vector<int>& colIndices,
                              Matrice* csrMatrix,
                              Matrice* cscMatrix,
                              int rowIndex,
                              int colIndex,
                              int& colPtrRef) {
    // colPtrRef = pointeur persistant pour cette colonne colIndex.
    // On démarre dessus, pas sur 0.

    int ix = 0;
    int iy = colPtrRef; // <-- RESTART ici
    int result = 0;

    if (rowIndices.empty() || colIndices.empty()) {
        return 0;
    }

    // merge classique CSR row vs CSC col,
    // mais la colonne repart à partir de colPtrRef
    while (ix < (int)rowIndices.size() && iy < (int)colIndices.size()) {
        int ridx = rowIndices[ix];
        int cidx = colIndices[iy];

        if (ridx == cidx) {
            // même index k => on multiplie A[i,k] * B[k,j]
            int rowValue = csrMatrix->GetPapier(rowIndex)->GetValeurById(ridx);
            int colValue = cscMatrix->GetPapier(cidx)->GetValeurById(colIndex);
            result += rowValue * colValue;

            ix++;
            iy++;
        } else if (ridx < cidx) {
            ix++;
        } else {
            // ridx > cidx
            iy++;
        }
    }

    // mise à jour du pointeur persistant :
    // On a avancé iy dans colIndices; on le sauvegarde
    colPtrRef = iy;

    return result;
}

Matrice* Matrice::MultiplyWithRestart(Matrice* m) {
    auto mtmp = new Matrice();
    mtmp->CopyAllPaper(this);

    const int nRows = this->GetNbPapier();
    const int nCols = m->GetNbPapier();

    // Préextraire toutes les lignes de A (CSR)
    std::vector<std::vector<int>> lines;     // pour chaque i: indices k tels que A[i,k] != 0
    lines.reserve(nRows);

    for (int i = 0; i < nRows; i++) {
        lines.push_back(this->GetPapier(i)->GetAllLignes());
        // (optionnel : tu as linesSize si tu veux, mais pas nécessaire pour l'algo)
    }

    // Préextraire toutes les colonnes de B (CSC)
    std::vector<std::vector<int>> columns;   // pour chaque j: indices k tels que B[k,j] != 0
    columns.reserve(nCols);

    for (int j = 0; j < nCols; j++) {
        columns.push_back(m->GetPapier(j)->GetAllColumns());
        // pareil : columnsSize si tu veux surveiller
    }

    // Voici le coeur du "restart":
    // col_ptrs[j] = où on reprend dans columns[j] la prochaine fois
    // On initialise tout à 0 (on n'a encore rien consommé des colonnes)
    std::vector<int> col_ptrs(nCols, 0);

    // On va parcourir les lignes i de this
    for (int i = 0; i < nRows; i++) {
        if (i % 10000 == 0) {
            std::cout << "Multiplication restart : "
                      << (i * 100) / std::max(1, nRows)
                      << " %" << std::endl;
        }

        // Raccourci vers les indices non nuls de la ligne i
        const std::vector<int>& rowIdx = lines[i];
        if (rowIdx.empty()) {
            continue;
        }

        // Pour chaque colonne j de m
        for (int j = 0; j < nCols; j++) {

            const std::vector<int>& colIdx = columns[j];
            if (colIdx.empty()) {
                continue;
            }

            // Pointeur persistant pour cette colonne j
            int& restartPtrForJ = col_ptrs[j];

            // Produit scalaire sparse entre:
            //   row i (rowIdx) et col j (colIdx),
            //   en démarrant l'exploration de la colonne à restartPtrForJ
            int value = dot_with_restart(
                rowIdx,
                colIdx,
                this,   // csrMatrix
                m,      // cscMatrix
                i,
                j,
                restartPtrForJ // <-- référence modifiable
            );

            if (value != 0) {
                mtmp->AddToLigne(value, i, j);
            }
        }
    }

    return mtmp;
}

Matrice* Matrice::MultiplyParallelWithRestart(Matrice* m){
    auto mtmp = new Matrice();
    mtmp->CopyAllPaper(this); // copie les métadonnées / papiers

    const int nRows = this->GetNbPapier();
    const int nCols = m->GetNbPapier();

    // Pré-extraire les lignes de A (CSR-like)
    std::vector<std::vector<int>> lines;
    lines.reserve(nRows);
    for (int i = 0; i < nRows; i++) {
        lines.push_back(this->GetPapier(i)->GetAllLignes());
    }

    // Pré-extraire les colonnes de B (CSC-like)
    std::vector<std::vector<int>> columns;
    columns.reserve(nCols);
    for (int j = 0; j < nCols; j++) {
        columns.push_back(m->GetPapier(j)->GetAllColumns());
    }

    // ------------------------------
    // PARALLÉLISATION
    // ------------------------------
    // On va paralléliser sur les lignes i.
    // Chaque thread :
    //   - a sa propre copie de col_ptrs
    //   - garde ses résultats locaux,
    // puis on fusionne après.
    // ------------------------------

    // Structure locale pour accumuler les valeurs non nulles calculées par un thread
    // clé: (i,j) encodé en 64 bits, valeur: somme entière
    auto encode_pair = [nCols](int i, int j) -> long long {
        return (static_cast<long long>(i) << 32) | static_cast<unsigned long long>(j);
    };

    // On va stocker les résultats de chaque thread dans un vecteur,
    // puis faire la fusion après la région parallèle.
    int max_threads = omp_get_max_threads();
    std::vector<std::map<long long,int>> thread_accums(max_threads);

    // On fait aussi une barre de progression grossière, mais attention :
    // impression depuis plusieurs threads = bruyant,
    // donc on ne l'affiche que depuis le thread master.
    // (optionnel, tu peux virer le bloc "if (omp_get_thread_num() == 0) { ... }")

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::map<long long,int>& local_accum = thread_accums[tid];

        // copie locale des pointeurs de redémarrage pour CHAQUE thread
        std::vector<int> col_ptrs_local(nCols, 0);

        // On parallélise sur i en for static
#pragma omp for schedule(static)
        for (int i = 0; i < nRows; i++) {

            // affichage de progression depuis le thread 0 uniquement
            if (tid == 0 && (i % 10000 == 0)) {
                std::cout << "Multiplication parallel-restart : "
                          << (i * 100) / std::max(1, nRows)
                          << " %" << std::endl;
            }

            const std::vector<int>& rowIdx = lines[i];
            if (rowIdx.empty()) {
                continue;
            }

            // Pour chaque colonne j
            for (int j = 0; j < nCols; j++) {
                const std::vector<int>& colIdx = columns[j];
                if (colIdx.empty()) {
                    continue;
                }

                int& restartPtrForJ = col_ptrs_local[j];

                // même appel que séquentiel
                int value = dot_with_restart(
                    rowIdx,
                    colIdx,
                    this,   // csrMatrix
                    m,      // cscMatrix
                    i,
                    j,
                    restartPtrForJ // pointeur modifié localement au thread
                );

                if (value != 0) {
                    long long key = encode_pair(i, j);
                    // additionne si (i,j) déjà vu par le même thread
                    local_accum[key] += value;
                }
            }
        } // end for i
    } // end parallel region

    // ------------------------------
    // Fusion des résultats dans mtmp
    // ------------------------------
    // Maintenant on insère toutes les contributions (i,j)->value finales.
    // On fusionne thread par thread. Ici on est séquentiel donc pas besoin de lock.
    for (int t = 0; t < max_threads; t++) {
        for (auto& kv : thread_accums[t]) {
            long long key = kv.first;
            int value = kv.second;
            if (value == 0) continue;

            int i = static_cast<int>(key >> 32);
            int j = static_cast<int>(key & 0xFFFFFFFF);

            mtmp->AddToLigne(value, i, j);
        }
    }

    return mtmp;
};

void Matrice::Add(Matrice *m) {
    for(int i = 0; i < m->GetNbPapier();i++){
        for(int j = 0; j< m->GetPapier(i)->GetNbLignes();j++){
            if(this->GetPapier(i)->IsInLignes(m->GetPapier(i)->GetLigne(j))){
                this->AddToLigne(m->GetPapier(i)->GetValeur(j),i,m->GetPapier(i)->GetLigne(j));
            }else{
                this->AddLigne(m->GetPapier(i)->GetValeur(j),i,m->GetPapier(i)->GetLigne(j));
                this->AddColumns(i,m->GetPapier(i)->GetLigne(j));
            }
        }
    }
    cout << "Addition effectuée" << endl;
}

Matrice::~Matrice() = default;

int Matrice::GetNbPapier() const {
    return this->nb_papier;
}

void Matrice::AddToLigne(int val, int idArticle, int idCite) {
    if(!this->GetPapier(idArticle)->AddToLigne(val,idCite)){
        this->AddColumns(idArticle,idCite);
    }
}

bool Matrice::IsInMatrice(int id) {
    return matrice.find(id) != matrice.end();
}

void Matrice::CopyAllPaper(Matrice *m) {
    for(auto it : m->matrice){
        this->AddPapier(it.first,it.second->GetTitre());
    }
}





