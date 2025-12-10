#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <iomanip>
#include <map>
#include <set>
#include <omp.h>

#include "Matrice.h"
#include "Papier.h"

using namespace std;

// ============================================================================
// STRUCTURES POUR RÉSULTATS EXPÉRIMENTAUX
// ============================================================================

struct ExperimentResult {
    string experiment_name;
    double beta_target;
    double beta_measured;
    double time_merge_ms;
    double time_restart_ms;
    double time_bitmap_only_ms;
    double time_combined_ms;
    double speedup_restart;
    double speedup_combined;
    double speedup_theory;
    double error_percent;
    int n, nnz;
    double rho_A, rho_B;
    int num_threads;
};

struct MatrixStats {
    int n, m, nnz;
    double density;
    double beta;
    double rho_A, rho_B;
    double rho_unique_B;
    int avg_row_nnz;
    int avg_col_nnz;
};

// ============================================================================
// CALCUL DE BETA (Definition IV.1 de l'article)
// ============================================================================

double compute_jaccard(const vector<int>& set1, const vector<int>& set2) {
    if (set1.empty() && set2.empty()) return 1.0;
    
    // Intersection
    vector<int> intersection;
    set_intersection(set1.begin(), set1.end(),
                    set2.begin(), set2.end(),
                    back_inserter(intersection));
    
    // Union
    vector<int> union_set;
    set_union(set1.begin(), set1.end(),
             set2.begin(), set2.end(),
             back_inserter(union_set));
    
    if (union_set.empty()) return 0.0;
    return static_cast<double>(intersection.size()) / union_set.size();
}

double compute_beta(Matrice* B) {
    int n = B->GetNbPapier();
    if (n <= 1) return 0.0;
    
    double sum_jaccard = 0.0;
    int valid_pairs = 0;
    
    // Pour chaque paire de colonnes consécutives (dans l'ordre des IDs)
    vector<int> sorted_ids;
    for (int j = 0; j < n; j++) {
        if (B->IsInMatrice(j)) {
            sorted_ids.push_back(j);
        }
    }
    sort(sorted_ids.begin(), sorted_ids.end());
    
    for (size_t i = 0; i < sorted_ids.size() - 1; i++) {
        int j = sorted_ids[i];
        int j_next = sorted_ids[i + 1];
        
        // Obtenir les indices non-zéros de chaque colonne
        vector<int> col_j = B->GetPapier(j)->GetAllColumns();
        vector<int> col_j_next = B->GetPapier(j_next)->GetAllColumns();
        
        sort(col_j.begin(), col_j.end());
        sort(col_j_next.begin(), col_j_next.end());
        
        double jaccard = compute_jaccard(col_j, col_j_next);
        sum_jaccard += jaccard;
        valid_pairs++;
    }
    
    return (valid_pairs > 0) ? (sum_jaccard / valid_pairs) : 0.0;
}

// ============================================================================
// CALCUL DE rho_unique_B (Definition IV.2)
// ============================================================================

int compute_rho_unique_B(Matrice* B) {
    set<int> unique_indices;
    int n = B->GetNbPapier();
    
    for (int j = 0; j < n; j++) {
        if (B->IsInMatrice(j)) {
            vector<int> col_indices = B->GetPapier(j)->GetAllColumns();
            unique_indices.insert(col_indices.begin(), col_indices.end());
        }
    }
    
    return unique_indices.size();
}

// ============================================================================
// STATISTIQUES DE MATRICE
// ============================================================================

MatrixStats compute_matrix_stats(Matrice* A, Matrice* B) {
    MatrixStats stats;
    stats.n = A->GetNbPapier();
    stats.m = A->GetNbPapier();
    stats.nnz = A->GetNnz();
    stats.density = static_cast<double>(stats.nnz) / (stats.n * stats.n);
    stats.beta = compute_beta(B);
    stats.rho_unique_B = compute_rho_unique_B(B);
    
    // Calculer densités moyennes
    int total_row_nnz = 0;
    int total_col_nnz = 0;
    int count = 0;
    
    for (int i = 0; i < stats.n; i++) {
        if (A->IsInMatrice(i)) {
            total_row_nnz += A->GetPapier(i)->GetNbLignes();
            total_col_nnz += A->GetPapier(i)->GetNbColumns();
            count++;
        }
    }
    
    stats.avg_row_nnz = (count > 0) ? (total_row_nnz / count) : 0;
    stats.avg_col_nnz = (count > 0) ? (total_col_nnz / count) : 0;
    stats.rho_A = stats.avg_row_nnz;
    stats.rho_B = stats.avg_col_nnz;
    
    return stats;
}

// ============================================================================
// GÉNÉRATEUR DE MATRICES SYNTHÉTIQUES AVEC β CONTRÔLÉ
// ============================================================================
// Crée une Matrice à partir d'une liste d'arêtes dirigées
// nodes[i] = nom du noeud i
// edges = liste de paires (src, dst) où src et dst sont des indices dans nodes
Matrice* build_custom_matrix(
    const vector<string>& nodes,
    const vector<pair<int,int>>& edges
) {
    Matrice* M = new Matrice();

    // 1. Déclarer tous les noeuds / colonnes
    for (int i = 0; i < (int)nodes.size(); i++) {
        M->AddPapier(i, nodes[i]);
    }

    // 2. Ajouter les arêtes comme des 1 dans la matrice clairsemée
    for (auto [src, dst] : edges) {
        // Mettre valeur 1 à la position (src, dst)
        M->AddLigne(1, src, dst);
        M->AddColumns(src, dst);
    }

    return M;
}

pair<vector<string>, vector<pair<int, int>>> load_snap_dataset(const string& filename) {
    vector<string> nodes;
    vector<pair<int, int>> edges;
    unordered_map<string, int> node_to_id;
    
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Impossible d'ouvrir le fichier: " + filename);
    }
    
    string line;
    int current_id = 0;
    
    while (getline(file, line)) {
        // Ignorer les lignes vides et les commentaires
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        istringstream iss(line);
        string node1, node2;
        
        // Lire les deux nœuds de l'arête
        if (!(iss >> node1 >> node2)) {
            continue; // Ligne mal formée, on l'ignore
        }
        
        // Ajouter node1 s'il n'existe pas
        if (node_to_id.find(node1) == node_to_id.end()) {
            node_to_id[node1] = current_id++;
            nodes.push_back(node1);
        }
        
        // Ajouter node2 s'il n'existe pas
        if (node_to_id.find(node2) == node_to_id.end()) {
            node_to_id[node2] = current_id++;
            nodes.push_back(node2);
        }
        
        // Ajouter l'arête avec les IDs numériques
        edges.push_back({node_to_id[node1], node_to_id[node2]});
    }
    
    file.close();
    return {nodes, edges};
}

Matrice* generate_synthetic_matrix(int n, double density, double beta_target, int seed = 42) {
    mt19937 rng(seed);
    
    Matrice* M = new Matrice();
    
    // Créer tous les papiers
    for (int i = 0; i < n; i++) {
        M->AddPapier(i, "Paper_" + to_string(i));
    }
    
    int nnz_per_col = static_cast<int>(n * density);
    if (nnz_per_col < 1) nnz_per_col = 1;
    
    // Colonne 0 : aléatoire
    vector<int> prev_col_indices;
    for (int k = 0; k < nnz_per_col; k++) {
        int row = uniform_int_distribution<>(0, n-1)(rng);
        prev_col_indices.push_back(row);
    }
    sort(prev_col_indices.begin(), prev_col_indices.end());
    prev_col_indices.erase(unique(prev_col_indices.begin(), prev_col_indices.end()), 
                           prev_col_indices.end());
    
    // Ajouter les entrées pour la colonne 0
    for (int row : prev_col_indices) {
        M->AddLigne(1, row, 0);
        M->AddColumns(row, 0);
    }
    
    // Colonnes suivantes : corrélées selon beta_target
    for (int j = 1; j < n; j++) {
        vector<int> curr_col_indices;
        
        if (prev_col_indices.empty()) {
            // Si colonne précédente vide, générer aléatoirement
            for (int k = 0; k < nnz_per_col; k++) {
                int row = uniform_int_distribution<>(0, n-1)(rng);
                curr_col_indices.push_back(row);
            }
        } else {
            // Garder beta_target * |prev| indices
            int keep_count = static_cast<int>(beta_target * prev_col_indices.size());
            shuffle(prev_col_indices.begin(), prev_col_indices.end(), rng);
            
            for (int k = 0; k < keep_count && k < prev_col_indices.size(); k++) {
                curr_col_indices.push_back(prev_col_indices[k]);
            }
            
            // Ajouter (1 - beta_target) * |prev| nouveaux indices
            int new_count = prev_col_indices.size() - keep_count;
            set<int> existing(curr_col_indices.begin(), curr_col_indices.end());
            
            while (curr_col_indices.size() < prev_col_indices.size() && 
                   curr_col_indices.size() < n) {
                int row = uniform_int_distribution<>(0, n-1)(rng);
                if (existing.find(row) == existing.end()) {
                    curr_col_indices.push_back(row);
                    existing.insert(row);
                }
            }
        }
        
        // Trier et supprimer doublons
        sort(curr_col_indices.begin(), curr_col_indices.end());
        curr_col_indices.erase(unique(curr_col_indices.begin(), curr_col_indices.end()), 
                              curr_col_indices.end());
        
        // Ajouter les entrées
        for (int row : curr_col_indices) {
            M->AddLigne(1, row, j);
            M->AddColumns(row, j);
        }
        
        prev_col_indices = curr_col_indices;
    }
    
    return M;
}

// ============================================================================
// FLUSH CACHE (pour benchmarks reproductibles)
// ============================================================================

void flush_cache() {
    const size_t size = 20 * 1024 * 1024; // 20 MB
    vector<char> dummy(size);
    for (size_t i = 0; i < size; i++) {
        dummy[i] = static_cast<char>(i);
    }
    // Force l'utilisation
    volatile char temp = dummy[size/2];
}

// ============================================================================
// BENCHMARKING FUNCTIONS
// ============================================================================

double benchmark_multiply_classic(Matrice* A, Matrice* B, int num_runs = 3) {
    vector<double> times;
    
    for (int run = 0; run < num_runs; run++) {
        flush_cache();
        
        auto start = chrono::high_resolution_clock::now();
        Matrice* C = A->Multiply(B);
        auto end = chrono::high_resolution_clock::now();

        C->ExportToJson("temp_classic.json");
        
        delete C;
        
        double ms = chrono::duration<double, milli>(end - start).count();
        times.push_back(ms);
    }
    
    // Retourner la médiane
    sort(times.begin(), times.end());
    return times[times.size() / 2];
}


double benchmark_multiply_restart(Matrice* A, Matrice* B, int num_runs = 3) {
    vector<double> times;
    
    for (int run = 0; run < num_runs; run++) {
        flush_cache();
        
        auto start = chrono::high_resolution_clock::now();
        Matrice* C = A->MultiplyWithRestart(B);
        auto end = chrono::high_resolution_clock::now();

        C->ExportToJson("temp_r.json");
        
        delete C;
        
        double ms = chrono::duration<double, milli>(end - start).count();
        times.push_back(ms);
    }
    
    sort(times.begin(), times.end());
    return times[times.size() / 2];
}

double benchmark_multiply(Matrice* A, Matrice* B, int num_runs = 3) {
    vector<double> times;
    
    for (int run = 0; run < num_runs; run++) {
        flush_cache();
        
        auto start = chrono::high_resolution_clock::now();
        Matrice* C = A->MultiplyNaive(B);
        auto end = chrono::high_resolution_clock::now();

        C->ExportToJson("temp_r.json");
        
        delete C;
        
        double ms = chrono::duration<double, milli>(end - start).count();
        times.push_back(ms);
    }
    
    sort(times.begin(), times.end());
    return times[times.size() / 2];
}

double benchmark_multiply_parallel_restart(Matrice* A, Matrice* B, int num_threads, int num_runs = 3) {
    vector<double> times;
    
    omp_set_num_threads(num_threads);
    
    for (int run = 0; run < num_runs; run++) {
        flush_cache();
        
        auto start = chrono::high_resolution_clock::now();
        Matrice* C = A->MultiplyParallelWithRestart(B);
        auto end = chrono::high_resolution_clock::now();
        
        delete C;
        
        double ms = chrono::duration<double, milli>(end - start).count();
        times.push_back(ms);
    }
    
    sort(times.begin(), times.end());
    return times[times.size() / 2];
}

// ============================================================================
// EXPÉRIENCE 1 : Q1 - VALIDATION β ⇒ SPEEDUP
// ============================================================================

void experiment_Q1_beta_vs_speedup(const string& output_file) {
    cout << "\n=== EXP-Q1: Beta vs Speedup ===" << endl;
    
    ofstream csv(output_file);
    csv << "beta_target,beta_measured,speedup_measured,speedup_theory,error_pct,time_merge_ms,time_restart_ms,n,nnz" << endl;
    
    int n = 4096;
    double density = 0.01;
    int num_trials = 10;
    
    //vector<double> beta_targets = {0.0, 0.1};
    //vector<double> beta_targets = {0.2, 0.3};
    //vector<double> beta_targets = {0.4, 0.5};
    //vector<double> beta_targets = {0.6, 0.7};
    //vector<double> beta_targets = {0.8, 0.9};
    //vector<double> beta_targets = {1.0};
    vector<double> beta_targets = {0.0,0.1,0.2, 0.3,0.4,0.5,0.6, 0.7,0.8, 0.9,0.99};

    for (double beta_target : beta_targets) {
        cout << "\nTesting beta = " << beta_target << endl;
        
        vector<double> speedups_measured;
        vector<double> speedups_theory;
        vector<double> betas_measured;
        vector<double> times_merge;
        vector<double> times_restart;
        
        for (int trial = 0; trial < num_trials; trial++) {
            // Générer matrices
            Matrice* A = generate_synthetic_matrix(n, density, beta_target, 42 + trial);
            Matrice* B = generate_synthetic_matrix(n, density, beta_target, 1000 + trial);
            
            // Mesurer beta réel
            double beta_actual = compute_beta(B);
            
            // Benchmark
            double t_merge = benchmark_multiply_classic(A, B, 3);
            double t_restart = benchmark_multiply_restart(A, B, 3);
            
            double speedup_measured = t_merge / t_restart;
            
            // Théorie (Theorem IV.3)
            MatrixStats stats = compute_matrix_stats(A, B);
            double rho_A = stats.rho_A;
            double rho_B = stats.rho_B;
            double speedup_theory = (rho_A + rho_B) / (rho_A + rho_B * (1.0 - beta_actual));
            
            betas_measured.push_back(beta_actual);
            speedups_measured.push_back(speedup_measured);
            speedups_theory.push_back(speedup_theory);
            times_merge.push_back(t_merge);
            times_restart.push_back(t_restart);
            
            delete A;
            delete B;
            
            cout << "  Trial " << (trial+1) << ": β_actual=" << beta_actual 
                 << ", speedup=" << speedup_measured << endl;
        }
        
        // Calculer moyennes
        double avg_beta = accumulate(betas_measured.begin(), betas_measured.end(), 0.0) / betas_measured.size();
        double avg_speedup_meas = accumulate(speedups_measured.begin(), speedups_measured.end(), 0.0) / speedups_measured.size();
        double avg_speedup_theory = accumulate(speedups_theory.begin(), speedups_theory.end(), 0.0) / speedups_theory.size();
        double avg_time_merge = accumulate(times_merge.begin(), times_merge.end(), 0.0) / times_merge.size();
        double avg_time_restart = accumulate(times_restart.begin(), times_restart.end(), 0.0) / times_restart.size();
        
        double error_pct = abs(avg_speedup_meas - avg_speedup_theory) / avg_speedup_theory * 100.0;
        
        csv << beta_target << "," << avg_beta << "," << avg_speedup_meas << "," 
            << avg_speedup_theory << "," << error_pct << "," << avg_time_merge << "," 
            << avg_time_restart << "," << n << "," << (int)(n * n * density) << endl;
        
        cout << "  AVERAGE: β=" << avg_beta << ", speedup=" << avg_speedup_meas 
             << " (theory=" << avg_speedup_theory << ", error=" << error_pct << "%)" << endl;
    }
    
    csv.close();
    cout << "\nResults saved to: " << output_file << endl;
}

void experiment_Q1_2_beta_vs_speedup(const string& output_file) {
    cout << "\n=== EXP-Q1_2: Density vs Speedup ===" << endl;
    
    ofstream csv(output_file);
    csv << "beta_target,beta_measured,speedup_measured,speedup_theory,error_pct,time_merge_ms,time_restart_ms,n,nnz,density" << endl;
    
    int n = 4096;
    vector<double> densities = {0.08, 0.09,0.099};
    double beta = 0.5;
    int num_trials = 3;
    
    for (double density : densities) {
        cout << "\nTesting density = " << density << endl;
        
        vector<double> speedups_measured;
        vector<double> speedups_theory;
        vector<double> betas_measured;
        vector<double> times_merge;
        vector<double> times_restart;
        
        for (int trial = 0; trial < num_trials; trial++) {
            // Générer matrices
            auto A = std::unique_ptr<Matrice>(generate_synthetic_matrix(n, density, beta, 42 + trial));
            auto B = std::unique_ptr<Matrice>(generate_synthetic_matrix(n, density, beta, 42 + trial));
            
            // Mesurer beta réel
            double beta_actual = compute_beta(B.get());
            
            // Benchmark
            double t_merge = benchmark_multiply_classic(A.get(), B.get(), 3);
            double t_restart = benchmark_multiply_restart(A.get(), B.get(), 3);
            
            double speedup_measured = t_merge / t_restart;
            
            // Théorie (Theorem IV.3)
            MatrixStats stats = compute_matrix_stats(A.get(), B.get());
            double rho_A = stats.rho_A;
            double rho_B = stats.rho_B;
            double speedup_theory = (rho_A + rho_B) / (rho_A + rho_B * (1.0 - beta_actual));
            
            betas_measured.push_back(beta_actual);
            speedups_measured.push_back(speedup_measured);
            speedups_theory.push_back(speedup_theory);
            times_merge.push_back(t_merge);
            times_restart.push_back(t_restart);

            
            cout << "  Trial " << (trial+1) << ": β_actual=" << beta_actual 
                 << ", speedup=" << speedup_measured << endl;
        }
        
        // Calculer moyennes
        double avg_beta = accumulate(betas_measured.begin(), betas_measured.end(), 0.0) / betas_measured.size();
        double avg_speedup_meas = accumulate(speedups_measured.begin(), speedups_measured.end(), 0.0) / speedups_measured.size();
        double avg_speedup_theory = accumulate(speedups_theory.begin(), speedups_theory.end(), 0.0) / speedups_theory.size();
        double avg_time_merge = accumulate(times_merge.begin(), times_merge.end(), 0.0) / times_merge.size();
        double avg_time_restart = accumulate(times_restart.begin(), times_restart.end(), 0.0) / times_restart.size();
        
        double error_pct = abs(avg_speedup_meas - avg_speedup_theory) / avg_speedup_theory * 100.0;
        
        csv << beta << "," << avg_beta << "," << avg_speedup_meas << "," 
            << avg_speedup_theory << "," << error_pct << "," << avg_time_merge << "," 
            << avg_time_restart << "," << n << "," << (int)(n * n * density) << "," << density << endl;
        
        cout << "  AVERAGE: β=" << avg_beta << ", speedup=" << avg_speedup_meas 
             << " (theory=" << avg_speedup_theory << ", error=" << error_pct << "%)" << endl;
    }
    
    csv.close();
    cout << "\nResults saved to: " << output_file << endl;
}

void experiment_Q6_beta_vs_speedup(const string& output_file) {
    cout << "\n=== EXP-Q6: Beta vs Speedup 3Dimensionnal===" << endl;
    
    ofstream csv(output_file);
    csv << "beta_target,beta_measured,speedup_measured,speedup_theory,error_pct,time_merge_ms,time_restart_ms,n,nnz,density" << endl;
    
    int num_trials = 3;
    

    //vector<int> n_targets = {32,64,128,256,512,1024,2048,4096};
    vector<int> n_targets = {4096};
    //vector<double> beta_targets = {0.0,0.1,0.2, 0.3,0.4,0.5,0.6, 0.7,0.8, 0.9,0.99};
    vector<double> beta_targets = {0.5,0.6,0.7,0.8, 0.9,0.99};
    vector<double> density_targets = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
    //vector<double> density_targets = {0.08,0.09,0.1};

    for (double beta_target : beta_targets) {
        for (double density_target : density_targets){
            for (int n_target : n_targets){
                cout << "\nTesting beta = " << beta_target << endl;
                
                vector<double> speedups_measured;
                vector<double> speedups_theory;
                vector<double> betas_measured;
                vector<double> times_merge;
                vector<double> times_restart;
                
                for (int trial = 0; trial < num_trials; trial++) {
                    // Générer matrices
                    auto A = std::unique_ptr<Matrice>(generate_synthetic_matrix(n_target, density_target, beta_target, 42 + trial));
                    auto B = std::unique_ptr<Matrice>(generate_synthetic_matrix(n_target, density_target, beta_target, 42 + trial));
                    
                    // Mesurer beta réel
                    double beta_actual = compute_beta(B.get());
                    
                    // Benchmark
                    double t_merge = benchmark_multiply_classic(A.get(), B.get(), 3);
                    double t_restart = benchmark_multiply_restart(A.get(), B.get(), 3);
                    
                    double speedup_measured = t_merge / t_restart;
                    
                    // Théorie (Theorem IV.3)
                    MatrixStats stats = compute_matrix_stats(A.get(), B.get());
                    double rho_A = stats.rho_A;
                    double rho_B = stats.rho_B;
                    double speedup_theory = (rho_A + rho_B) / (rho_A + rho_B * (1.0 - beta_actual));
                    
                    betas_measured.push_back(beta_actual);
                    speedups_measured.push_back(speedup_measured);
                    speedups_theory.push_back(speedup_theory);
                    times_merge.push_back(t_merge);
                    times_restart.push_back(t_restart);
                
                    
                    cout << "  Trial " << (trial+1) << ": β_actual=" << beta_actual 
                        << ", speedup=" << speedup_measured << ", density= "
                        << density_target <<", n=" << n_target << endl;
                }
                
                // Calculer moyennes
                double avg_beta = accumulate(betas_measured.begin(), betas_measured.end(), 0.0) / betas_measured.size();
                double avg_speedup_meas = accumulate(speedups_measured.begin(), speedups_measured.end(), 0.0) / speedups_measured.size();
                double avg_speedup_theory = accumulate(speedups_theory.begin(), speedups_theory.end(), 0.0) / speedups_theory.size();
                double avg_time_merge = accumulate(times_merge.begin(), times_merge.end(), 0.0) / times_merge.size();
                double avg_time_restart = accumulate(times_restart.begin(), times_restart.end(), 0.0) / times_restart.size();
                
                double error_pct = abs(avg_speedup_meas - avg_speedup_theory) / avg_speedup_theory * 100.0;
                
                csv << beta_target << "," << avg_beta << "," << avg_speedup_meas << "," 
                    << avg_speedup_theory << "," << error_pct << "," << avg_time_merge << "," 
                    << avg_time_restart << "," << n_target << "," << (int)(n_target * n_target * density_target) << "," << density_target<<  endl;
                
                cout << "  AVERAGE: β=" << avg_beta << ", speedup=" << avg_speedup_meas 
                    << " (theory=" << avg_speedup_theory << ", error=" << error_pct << "%)" << endl;

            }

        }
        
    }
    
    csv.close();
    cout << "\nResults saved to: " << output_file << endl;
}

// ============================================================================
// EXPÉRIENCE 2 : Q2 - ABLATION STUDY
// ============================================================================

void experiment_Q2_ablation_study(const string& output_file) {
    cout << "\n=== EXP-Q2: Ablation Study ===" << endl;
    
    ofstream csv(output_file);
    csv << "matrix_type,beta,speedup_restart_only,speedup_combined,composition_ratio,time_merge_ms,time_restart_ms,time_combined_ms" << endl;
    csv << "matrix_type,beta,time_naive_ms" << endl;
    
    struct TestMatrix {
        string name;
        string filename;
    };
    
    vector<TestMatrix> test_cases = {
        {"Enron", "/home/ciad/project_anstett/BENCH_NZIM_SPEGMM/build/email-Enron.txt"},
        {"Stan", "/home/ciad/project_anstett/BENCH_NZIM_SPEGMM/build/web-Stanford.txt"},
        {"Wiki", "/home/ciad/project_anstett/BENCH_NZIM_SPEGMM/build/wiki-Talk.txt"},
        {"ca-GrQc", "/home/ciad/project_anstett/BENCH_NZIM_SPEGMM/build/ca-GrQc.txt"},
        {"wiki-Vote", "/home/ciad/project_anstett/BENCH_NZIM_SPEGMM/build/wiki-Vote.txt"},
        {"p2p-Gnutella08", "/home/ciad/project_anstett/BENCH_NZIM_SPEGMM/build/p2p-Gnutella08.txt"},
        {"facebook_combined", "/home/ciad/project_anstett/BENCH_NZIM_SPEGMM/build/facebook_combined.txt"}
    };
    
    int n = 8;
    double density = 0.01;
    
    for (const auto& test : test_cases) {
        
        auto [nodes, edges] = load_snap_dataset(test.filename);
        
        Matrice* A = build_custom_matrix(nodes, edges);
        Matrice* B = build_custom_matrix(nodes, edges);
        
        double beta_measured  = compute_beta(B);

        cout << "\nTesting " << test.name << " (beta=" << beta_measured << ")" << endl;
        
        // Benchmark variantes
        double t_restart =  benchmark_multiply_classic(A, B, 2);
        double t_combined = benchmark_multiply_restart(A, B, 2);
        double t_naif = benchmark_multiply(A, B, 5);
        
        double composition =   t_restart / t_combined;
        
        csv << test.name << "," << beta_measured << "," 
            << "," << composition << "," 
            << t_restart << "," << t_combined << "," << t_naif << endl;
        
        cout << "  Restart: " << t_restart << " ms " << endl;
        cout << "  Combined: " << t_combined << " ms (speedup=" << composition << "x)" << endl;
        
        delete A;
        delete B;
    }
    
    csv.close();
    cout << "\nResults saved to: " << output_file << endl;
}

// ============================================================================
// EXPÉRIENCE 3 : Q3 - AMORTIZATION
// ============================================================================

void experiment_Q3_amortization(const string& output_file) {
    cout << "\n=== EXP-Q3: Bitmap Construction Amortization ===" << endl;
    
    ofstream csv(output_file);
    csv << "k_iterations,time_per_iter_merge_ms,time_per_iter_restart_ms,speedup,breakeven" << endl;
    
    int n = 4096;
    double density = 0.01;
    double beta = 0.6;
    
    Matrice* A = generate_synthetic_matrix(n, density, beta, 42);
    Matrice* B = generate_synthetic_matrix(n, density, beta, 1000);
    
    vector<int> k_values = {1, 2, 3, 5, 10, 15, 20};
    
    for (int k : k_values) {
        cout << "\nTesting k=" << k << " iterations" << endl;
        
        // Simuler k multiplications successives
        double total_time_merge = 0.0;
        double total_time_restart = 0.0;
        
        for (int iter = 0; iter < k; iter++) {
            total_time_merge += benchmark_multiply_classic(A, B, 1);
            total_time_restart += benchmark_multiply_restart(A, B, 1);
        }
        
        double time_per_iter_merge = total_time_merge / k;
        double time_per_iter_restart = total_time_restart / k;
        double speedup = time_per_iter_merge / time_per_iter_restart;
        bool breakeven = (speedup >= 1.0);
        
        csv << k << "," << time_per_iter_merge << "," << time_per_iter_restart << "," 
            << speedup << "," << (breakeven ? "1" : "0") << endl;
        
        cout << "  Per-iteration time: Merge=" << time_per_iter_merge << " ms, "
             << "Restart=" << time_per_iter_restart << " ms (speedup=" << speedup << "x)" << endl;
        cout << "  Breakeven: " << (breakeven ? "YES" : "NO") << endl;
    }
    
    delete A;
    delete B;
    
    csv.close();
    cout << "\nResults saved to: " << output_file << endl;
}

// ============================================================================
// EXPÉRIENCE 4 : Q4 - SCALABILITY
// ============================================================================

void experiment_Q4_scalability(const string& output_file) {
    cout << "\n=== EXP-Q4: Scalability Analysis ===" << endl;
    
    ofstream csv(output_file);
    csv << "num_threads,time_ms,speedup,efficiency,matrix_size" << endl;
    
    int n = 8192;
    double density = 0.005;
    double beta = 0.6;
    
    
    vector<int> thread_counts = {1, 2, 4, 8, 16};
    //For weak scaling
    vector<int> actual_thread_matrix_sizes = {2048, 2896,4096,5776,8192};
    size_t index = 0;
    double t_sequential = 0.0;
    
    for (int num_threads : thread_counts) {
        cout << "\nTesting with " << num_threads << " threads" << endl;

        Matrice* A = generate_synthetic_matrix(actual_thread_matrix_sizes[index], density, beta, 42);
        Matrice* B = generate_synthetic_matrix(actual_thread_matrix_sizes[index], density, beta, 1000);

        
        double t_parallel = benchmark_multiply_parallel_restart(A, B, num_threads, 5);
        
        if (num_threads == 1) {
            t_sequential = t_parallel;
        }
        
        double speedup = t_sequential / t_parallel;
        double efficiency = speedup / num_threads;
        
        csv << num_threads << "," << t_parallel << "," << speedup << "," 
            << efficiency << "," << actual_thread_matrix_sizes[index] << endl;

        index++;
        
        cout << "  Time: " << t_parallel << " ms" << endl;
        cout << "  Speedup: " << speedup << "x (efficiency: " << (efficiency*100) << "%)" << endl;
        delete A;
        delete B;
    }
    
    
    
    csv.close();
    cout << "\nResults saved to: " << output_file << endl;
}

// ============================================================================
// EXPÉRIENCE 5 : Q5 - UNFAVORABLE CASES
// ============================================================================

void experiment_Q5_unfavorable_cases(const string& output_file) {
    cout << "\n=== EXP-Q5: Unfavorable Cases ===" << endl;
    
    ofstream csv(output_file);
    csv << "case_name,n,beta,speedup,time_merge_ms,time_restart_ms,favorable" << endl;
    
    // Cas 1: Matrice aléatoire (β ≈ 0)
    cout << "\nCase 1: Random matrix" << endl;
    {
        int n = 4096;
        Matrice* A = generate_synthetic_matrix(n, 0.01, 0.0, 42);
        Matrice* B = generate_synthetic_matrix(n, 0.01, 0.0, 1000);
        
        double beta = compute_beta(B);
        double t_merge = benchmark_multiply_classic(A, B, 5);
        double t_restart = benchmark_multiply_restart(A, B, 5);
        double speedup = t_merge / t_restart;
        
        csv << "Random," << n << "," << beta << "," << speedup << "," 
            << t_merge << "," << t_restart << "," << (speedup > 1.05 ? "1" : "0") << endl;
        
        cout << "  β=" << beta << ", speedup=" << speedup << endl;
        
        delete A;
        delete B;
    }
    
    // Cas 2: Petites matrices
    cout << "\nCase 2: Small matrices" << endl;
    vector<int> sizes = {64, 128, 256, 512, 1024, 2048};
    for (int n : sizes) {
        Matrice* A = generate_synthetic_matrix(n, 0.01, 0.6, 42);
        Matrice* B = generate_synthetic_matrix(n, 0.01, 0.6, 1000);
        
        double beta = compute_beta(B);
        double t_merge = benchmark_multiply_classic(A, B, 3);
        double t_restart = benchmark_multiply_restart(A, B, 3);
        double speedup = t_merge / t_restart;
        
        csv << "Small_n" << n << "," << n << "," << beta << "," << speedup << "," 
            << t_merge << "," << t_restart << "," << (speedup > 1.05 ? "1" : "0") << endl;
        
        cout << "  n=" << n << ", β=" << beta << ", speedup=" << speedup << endl;
        
        delete A;
        delete B;
    }
    
    // Cas 3: Beta modéré
    cout << "\nCase 3: Moderate beta" << endl;
    {
        int n = 4096;
        Matrice* A = generate_synthetic_matrix(n, 0.01, 0.3, 42);
        Matrice* B = generate_synthetic_matrix(n, 0.01, 0.3, 1000);
        
        double beta = compute_beta(B);
        double t_merge = benchmark_multiply_classic(A, B, 5);
        double t_restart = benchmark_multiply_restart(A, B, 5);
        double speedup = t_merge / t_restart;
        
        csv << "Moderate_beta," << n << "," << beta << "," << speedup << "," 
            << t_merge << "," << t_restart << "," << (speedup > 1.05 ? "1" : "0") << endl;
        
        cout << "  β=" << beta << ", speedup=" << speedup << endl;
        
        delete A;
        delete B;
    }
    
    csv.close();
    cout << "\nResults saved to: " << output_file << endl;
}

// ============================================================================
// FONCTION PRINCIPALE
// ============================================================================

int main(int argc, char* argv[]) {
    cout << "========================================" << endl;
    cout << "NZIM SpGEMM Experimental Framework" << endl;
    cout << "Article: Restart-Based Sparse Matrix Multiplication" << endl;
    cout << "========================================" << endl;
    
    string output_dir = "/home/ciad/project_anstett/BENCH_NZIM_SPEGMM/result/";
    
    // Lancer toutes les expériences
    cout << "\n>>> Starting experimental validation suite <<<\n" << endl;
    
    try {
        // Q1: Beta vs Speedup (Figure 2)
        //experiment_Q1_beta_vs_speedup(output_dir + "exp_q1_beta_speedup.csv");

        //Q1: Density vs Speedup (Figure 2.1)
        //experiment_Q1_2_beta_vs_speedup(output_dir + "exp_q1_2_density_speedup.csv");
        
        // Q2: Ablation Study (Table III)
        experiment_Q2_ablation_study(output_dir + "exp_q2_ablation.csv");
        
        // Q3: Amortization (Figure 3)
        //experiment_Q3_amortization(output_dir + "exp_q3_amortization.csv");
        
        // Q4: Scalability (Figure 4)
        //experiment_Q4_scalability(output_dir + "exp_q4_scalability.csv");
        
        // Q5: Unfavorable Cases (Table in Section VI)
        //experiment_Q5_unfavorable_cases(output_dir + "exp_q5_unfavorable.csv");

        // Q1: Beta vs Speedup (Figure 2)
        //experiment_Q6_beta_vs_speedup(output_dir + "exp_q6_beta_speedup.csv");
        
        cout << "\n========================================" << endl;
        cout << "ALL EXPERIMENTS COMPLETED SUCCESSFULLY" << endl;
        cout << "Results saved to: " << output_dir << endl;
        cout << "========================================" << endl;
        
    } catch (const exception& e) {
        cerr << "\nERROR: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
