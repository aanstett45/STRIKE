//
// Created by alex on 22/03/24.
//

#ifndef INA_NK_COR_PAPIER_H
#define INA_NK_COR_PAPIER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <roaring/roaring.hh>

using namespace std;

class Papier {
private :
    string title;
    int nnz;
    int nb_lignes;
    int nb_columns;
    RoaringBitmap Lignes;
    RoaringBitmap Columns;
    unordered_map<int,int> Valeurs;
public :
    explicit Papier(const string& titre);
    explicit  Papier(Papier * papier);
    void SetTitre(const string& titre);
    string GetTitre();
    int GetNbColumns() const;
    int GetNbLignes() const;
    int GetColumns(int id);
    int GetLigne(int id);
    int GetValeur(int id);
    const unordered_map<int,int>& getValeurs() const { return Valeurs; }
    int GetNnz() const;
    vector<int> GetAllValues();
    int GetValeurById(int id);

    void ensureLignesSize(int index);

    void ensureColumnsSize(int index);

    void AddLigne(int val, int id);
    void AddColumn(int id);
    bool AddToLigne(int val,int id);
    void DeleteLigne(int id);
    void DeleteColumn(int id);
    bool IsInLignes(int id);
    vector<int> GetAllLignes();
    vector<int> GetAllColumns();

    const RoaringBitmap &getLignesBitmap() const;

    const RoaringBitmap &getColumnsBitmap() const;

    ~Papier();
};


#endif //INA_NK_COR_PAPIER_H
