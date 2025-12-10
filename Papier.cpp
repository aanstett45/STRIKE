//
// Created by alex on 22/03/24.
//

#include <algorithm>
#include "Papier.h"

Papier::Papier(const string& titre): Lignes(), Columns() {
    this->title = titre;
    this->nb_lignes = 0;
    this->nb_columns = 0;
    this->nnz = 0;
}

Papier::Papier(Papier *papier) {
    title = papier->title;
    nb_lignes = papier->nb_lignes;
    nb_columns = papier->nb_columns;
    this->nnz = papier->nnz;

    Columns = papier->Columns;
    Lignes = papier->Lignes;
    Valeurs = papier->Valeurs;
}

void Papier::SetTitre(const string& titre) {
    this->title = titre;
}

int Papier::GetNnz() const {
    return nnz;
}

string Papier::GetTitre() {
    return this->title;
}

int Papier::GetColumns(int id) {
    int count = 0;
    for (auto i : Columns) {  // Roaring bitmaps can be iterated directly
        if (count == id)
            return static_cast<int>(i);
        count++;
    }
    return -1;
}

int Papier::GetLigne(int id) {
    int count = 0;
    for (auto i : Lignes) {
        if (count == id)
            return static_cast<int>(i);
        count++;
    }
    return -1;
}

int Papier::GetValeur(int id) {
    int actualIndex = GetLigne(id);
    return Valeurs.at(actualIndex);
}

vector<int> Papier::GetAllValues() {
    vector<int> ret;
    for (const auto& kv : Valeurs) {
        ret.push_back(kv.second);
    }
    return ret;
}

// Modifier GetValeurById pour unordered_map
int Papier::GetValeurById(int id) {
    auto it = Valeurs.find(id);
    if (it != Valeurs.end()) {
        return it->second;
    }
    return 0;  // Ou throw une exception selon ton besoin
}

void Papier::AddLigne(int val, int id) {
    bool alreadySet = Lignes.contains(id);
    Lignes.add(id);  // Add the bit to the bitmap
    if (!alreadySet) {
        nb_lignes++;
    }
    Valeurs[id] = val;
    nnz ++;
}

void Papier::AddColumn(int id) {
    bool alreadySet = Columns.contains(id);
    Columns.add(id);
    if (!alreadySet) {
        nb_columns++;
    }
}


bool Papier::AddToLigne(int val, int id) {
    if(IsInLignes(id)){
        Valeurs[id] += val;
        return true;
    } else {
        AddLigne(val, id);
        return false;
    }
}

void Papier::DeleteLigne(int id) {
    if(Lignes.contains(id)) {  // Check if id is in the bitmap
        Lignes.remove(id);     // Remove the id from the bitmap
        Valeurs.erase(id);
        nb_lignes--;
        nnz--;
    }
}

void Papier::DeleteColumn(int id) {
    if(Columns.contains(id)) {  // Check if id exists in the bitmap
        Columns.remove(id);     // Remove id from the bitmap
        nb_columns--;
    }
}

Papier::~Papier() = default;

int Papier::GetNbColumns() const {
    return nb_columns;
}

int Papier::GetNbLignes() const {
    return nb_lignes;
}

bool Papier::IsInLignes(int id) {
    return Lignes.contains(id);
}

vector<int> Papier::GetAllLignes() {
    vector<int> result;
    for (auto i : Lignes) {  // Roaring provides direct iteration over set bits
        result.push_back(static_cast<int>(i));
    }
    return result;
}

vector<int> Papier::GetAllColumns() {
    vector<int> result;
    for (auto i : Columns) {  // Roaring provides direct iteration over set bits
        result.push_back(static_cast<int>(i));
    }
    return result;
}

const RoaringBitmap& Papier::getLignesBitmap() const {
    return Lignes;
}

const RoaringBitmap& Papier::getColumnsBitmap() const {
    return Columns;
}



