#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>

#include "Matrice.h"
#include "Papier.h"

using namespace std;

// ============================================================================
// Calcul de β (copié du framework expérimental)
// ============================================================================

double compute_jaccard(const vector<int>& set1, const vector<int>& set2) {
    if (set1.empty() && set2.empty()) return 1.0;
    
    vector<int> intersection;
    set_intersection(set1.begin(), set1.end(),
                    set2.begin(), set2.end(),
                    back_inserter(intersection));
    
    vector<int> union_set;
    set_union(set1.begin(), set1.end(),
             set2.begin(), set2.end(),
             back_inserter(union_set));
    
    if (union_set.empty()) return 0.0;
    return static_cast<double>(intersection.size()) / union_set.size();
}

double compute_beta(Matrice* M) {
    int n = M->GetNbPapier();
    if (n <= 1) {
        cout << "Warning: Matrix has ≤ 1 columns, β undefined" << endl;
        return 0.0;
    }
    
    double sum_jaccard = 0.0;
    int valid_pairs = 0;
    
    vector<int> sorted_ids;
    for (int j = 0; j < n; j++) {
        if (M->IsInMatrice(j)) {
            sorted_ids.push_back(j);
        }
    }
    sort(sorted_ids.begin(), sorted_ids.end());
    
    cout << "Computing β over " << (sorted_ids.size() - 1) << " consecutive column pairs..." << endl;
    
    for (size_t i = 0; i < sorted_ids.size() - 1; i++) {
        int j = sorted_ids[i];
        int j_next = sorted_ids[i + 1];
        
        vector<int> col_j = M->GetPapier(j)->GetAllColumns();
        vector<int> col_j_next = M->GetPapier(j_next)->GetAllColumns();
        
        sort(col_j.begin(), col_j.end());
        sort(col_j_next.begin(), col_j_next.end());
        
        double jaccard = compute_jaccard(col_j, col_j_next);
        sum_jaccard += jaccard;
        valid_pairs++;
        
        // Afficher progression tous les 100 pairs
        if ((i + 1) % 100 == 0) {
            cout << "  Processed " << (i + 1) << " / " << (sorted_ids.size() - 1) << " pairs..." << endl;
        }
    }
    
    return (valid_pairs > 0) ? (sum_jaccard / valid_pairs) : 0.0;
}

// ============================================================================
// Statistiques détaillées
// ============================================================================

void print_detailed_stats(Matrice* M) {
    int n = M->GetNbPapier();
    int nnz = M->GetNnz();
    
    cout << "\n========================================" << endl;
    cout << "MATRIX STATISTICS" << endl;
    cout << "========================================" << endl;
    cout << "Dimensions: " << n << " × " << n << endl;
    cout << "Non-zeros: " << nnz << endl;
    cout << "Density: " << scientific << setprecision(6) 
         << (static_cast<double>(nnz) / (static_cast<long long>(n) * n)) << endl;
    
    // Statistiques par ligne/colonne
    vector<int> row_nnzs, col_nnzs;
    for (int i = 0; i < n; i++) {
        if (M->IsInMatrice(i)) {
            row_nnzs.push_back(M->GetPapier(i)->GetNbLignes());
            col_nnzs.push_back(M->GetPapier(i)->GetNbColumns());
        }
    }
    
    if (!row_nnzs.empty()) {
        sort(row_nnzs.begin(), row_nnzs.end());
        sort(col_nnzs.begin(), col_nnzs.end());
        
        double avg_row = accumulate(row_nnzs.begin(), row_nnzs.end(), 0.0) / row_nnzs.size();
        double avg_col = accumulate(col_nnzs.begin(), col_nnzs.end(), 0.0) / col_nnzs.size();
        
        cout << "\nRow statistics:" << endl;
        cout << "  Average nnz: " << avg_row << endl;
        cout << "  Median nnz: " << row_nnzs[row_nnzs.size() / 2] << endl;
        cout << "  Max nnz: " << row_nnzs.back() << endl;
        
        cout << "\nColumn statistics:" << endl;
        cout << "  Average nnz: " << avg_col << endl;
        cout << "  Median nnz: " << col_nnzs[col_nnzs.size() / 2] << endl;
        cout << "  Max nnz: " << col_nnzs.back() << endl;
    }
    
    cout << "========================================" << endl;
}

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char* argv[]) {
    cout << "========================================" << endl;
    cout << "BETA CALCULATOR TOOL" << endl;
    cout << "Computes column correlation β (Definition IV.1)" << endl;
    cout << "========================================\n" << endl;
    
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <matrix_file> [output_csv]" << endl;
        cerr << "\nSupported formats:" << endl;
        cerr << "  - JSON (DBLP format)" << endl;
        cerr << "  - CSV (id,references)" << endl;
        return 1;
    }
    
    string input_file = argv[1];
    string output_file = (argc >= 3) ? argv[2] : "";
    
    cout << "Loading matrix from: " << input_file << endl;
    
    Matrice* M = nullptr;
    
    try {
        // Détecter le format
        if (input_file.find(".json") != string::npos) {
            M = new Matrice(input_file);
            cout << "✓ Loaded JSON format" << endl;
        } else if (input_file.find(".csv") != string::npos) {
            // Nécessite un fichier dictionnaire
            string dict_file = input_file;
            dict_file.replace(dict_file.find(".csv"), 4, "_dict.csv");
            M = new Matrice(input_file, dict_file, 1024);
            cout << "✓ Loaded CSV format" << endl;
        } else {
            cerr << "Error: Unsupported file format" << endl;
            return 1;
        }
        
        // Afficher statistiques
        print_detailed_stats(M);
        
        // Calculer β
        cout << "\nComputing β (column correlation)..." << endl;
        auto start = chrono::high_resolution_clock::now();
        double beta = compute_beta(M);
        auto end = chrono::high_resolution_clock::now();
        double elapsed_ms = chrono::duration<double, milli>(end - start).count();
        
        cout << "\n========================================" << endl;
        cout << "RESULT" << endl;
        cout << "========================================" << endl;
        cout << "β (column correlation): " << fixed << setprecision(4) << beta << endl;
        cout << "Computation time: " << elapsed_ms << " ms" << endl;
        
        // Interpréter la valeur
        cout << "\nInterpretation:" << endl;
        if (beta < 0.15) {
            cout << "  → β < 0.15: NO BENEFIT expected from restart-based SpGEMM" << endl;
            cout << "  → Recommendation: Use merge-based or MKL" << endl;
        } else if (beta < 0.3) {
            cout << "  → 0.15 ≤ β < 0.3: MARGINAL BENEFIT (10-15% speedup)" << endl;
            cout << "  → Recommendation: Consider column reordering first" << endl;
        } else if (beta < 0.5) {
            cout << "  → 0.3 ≤ β < 0.5: MODERATE BENEFIT (15-25% speedup)" << endl;
            cout << "  → Recommendation: Restart-based can be beneficial" << endl;
        } else {
            cout << "  → β ≥ 0.5: SUBSTANTIAL BENEFIT (>30% speedup)" << endl;
            cout << "  → Recommendation: USE restart-based SpGEMM" << endl;
        }
        cout << "========================================" << endl;
        
        // Sauvegarder dans CSV si demandé
        if (!output_file.empty()) {
            ofstream csv(output_file);
            csv << "matrix_file,n,nnz,density,beta,computation_time_ms" << endl;
            csv << input_file << "," << M->GetNbPapier() << "," << M->GetNnz() << ","
                << scientific << (static_cast<double>(M->GetNnz()) / (static_cast<long long>(M->GetNbPapier()) * M->GetNbPapier())) << ","
                << fixed << beta << "," << elapsed_ms << endl;
            csv.close();
            cout << "\n✓ Results saved to: " << output_file << endl;
        }
        
        delete M;
        
    } catch (const exception& e) {
        cerr << "\nERROR: " << e.what() << endl;
        if (M) delete M;
        return 1;
    }
    
    return 0;
}
