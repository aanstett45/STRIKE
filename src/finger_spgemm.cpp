// ============================================================================
// finger_spgemm.cpp -- Persistent-finger sparse intersections
//
// Instrument expérimental pour le papier révisé. Autonome : C++17 + OpenMP,
// aucune dépendance externe (Roaring supprimé : les merges n'en ont pas
// besoin, et les valeurs sont désormais dans de vrais tableaux CSR/CSC --
// plus aucun unordered_map dans les noyaux).
//
// Contributions instrumentées :
//   A. Noyau d'intersection à finger persistant, GARDÉ (Théorème 2) :
//      - garde O(1) + rembobinage par recherche binaire (exact partout)
//      - avancée par galloping (coût ~ O(log saut + |intersection|))
//   B. SpGEMM masqué / comptage de triangles (C = A^2 .* A) :
//      - le masque supprime le terme m*p ; sans persistance la colonne j
//        est scannée deg(j) fois, avec persistance <= 1 fois au total
//      - baselines : Gustavson masqué, edge-merge, edge-gallop,
//        compact-forward (merge/gallop/finger)
//   C. Taux de rembobinage = prédicteur en ligne ; ordonnancements
//      naturel / RCM / degré ; sélecteur adaptatif par échantillonnage
//   D. SpGEMM complet (positionnement honnête) : Gustavson SPA vs merge
//      stateless vs merge gardé(+gallop), balayage beta + bandé,
//      vérification de correction systématique
//
// Ordonnancement OpenMP : les boucles paralleles utilisent schedule(runtime);
// fixer OMP_SCHEDULE=static (defaut recommande), dynamic,64 ou guided pour
// la comparaison d'ordonnancements demandee par les relecteurs. Chaque
// benchmark comporte desormais un echauffement non mesure et imprime sa
// dispersion (min/med/max) sur stdout.
// Compilation :  g++ -O3 -march=native -fopenmp -std=c++17 \
//                    finger_spgemm.cpp -o finger_spgemm
//
// Commandes :
//   finger_spgemm verify
//   finger_spgemm fullspgemm <out.csv> [n=2048] [trials=3] [runs=3]
//   finger_spgemm triangles  <out.csv> [--runs R] [--threads T] <graphe...>
//   finger_spgemm adaptive   <out.csv> [--sample-pct P] <graphe...>
//
// <graphe> : fichier.txt (liste d'arêtes SNAP, symétrisée automatiquement)
//            fichier.mtx (MatrixMarket coordinate, pattern ou valeurs)
//            gen:er:n:deg      (Erdos-Renyi, degré moyen deg)
//            gen:ba:n:m        (Barabasi-Albert, m arêtes/noeud, degrés biaisés)
//            gen:band:n:w      (graphe bande : i ~ i+1..i+w, régime favorable)
// ============================================================================

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <omp.h>

using namespace std;

static long long g_sink = 0; // anti-DCE pour les benchmarks

// ============================================================================
// STRUCTURES CSR / CSC
// ============================================================================

struct CSR {
    int n = 0, m = 0;            // lignes, colonnes
    vector<long long> indptr;    // taille n+1
    vector<int> idx;             // indices de colonne, triés par ligne
    vector<double> val;          // valeurs alignées sur idx
    long long nnz() const { return (long long)idx.size(); }
};

// CSC(B) est stocké comme CSR(B^T) : indptr par colonne, idx = indices de
// ligne triés, val alignées.
static CSR transpose(const CSR& A) {
    CSR T;
    T.n = A.m; T.m = A.n;
    T.indptr.assign((size_t)T.n + 1, 0);
    for (long long e = 0; e < A.nnz(); e++) T.indptr[(size_t)A.idx[e] + 1]++;
    for (int j = 0; j < T.n; j++) T.indptr[(size_t)j + 1] += T.indptr[j];
    T.idx.resize(A.nnz());
    T.val.resize(A.nnz());
    vector<long long> pos(T.indptr.begin(), T.indptr.end() - 1);
    for (int i = 0; i < A.n; i++) {
        for (long long e = A.indptr[i]; e < A.indptr[(size_t)i + 1]; e++) {
            long long p = pos[A.idx[e]]++;
            T.idx[p] = i;                 // i croissant => colonnes triées
            T.val[p] = A.val[e];
        }
    }
    return T;
}

// Construction depuis triplets (i,j,v) : tri, fusion des doublons (somme).
static CSR from_triplets(int n, int m, vector<tuple<int,int,double>>& t) {
    sort(t.begin(), t.end(),
         [](const auto& a, const auto& b) {
             return make_pair(get<0>(a), get<1>(a)) <
                    make_pair(get<0>(b), get<1>(b));
         });
    CSR A; A.n = n; A.m = m;
    A.indptr.assign((size_t)n + 1, 0);
    size_t k = 0;
    while (k < t.size()) {
        int i = get<0>(t[k]), j = get<1>(t[k]);
        double s = 0;
        while (k < t.size() && get<0>(t[k]) == i && get<1>(t[k]) == j)
            s += get<2>(t[k++]);
        if (s != 0.0) {
            A.idx.push_back(j);
            A.val.push_back(s);
            A.indptr[(size_t)i + 1]++;
        }
    }
    for (int i = 0; i < n; i++) A.indptr[(size_t)i + 1] += A.indptr[i];
    return A;
}

static bool same_csr(const CSR& A, const CSR& B, long long* nb_diff = nullptr,
                     double eps = 1e-9) {
    long long diff = 0;
    if (A.n != B.n || A.m != B.m) { if (nb_diff) *nb_diff = -1; return false; }
    for (int i = 0; i < A.n; i++) {
        long long a = A.indptr[i], ae = A.indptr[(size_t)i + 1];
        long long b = B.indptr[i], be = B.indptr[(size_t)i + 1];
        while (a < ae && b < be) {
            if (A.idx[a] == B.idx[b]) {
                if (fabs(A.val[a] - B.val[b]) > eps) diff++;
                a++; b++;
            } else if (A.idx[a] < B.idx[b]) { diff++; a++; }
            else { diff++; b++; }
        }
        diff += (ae - a) + (be - b);
    }
    if (nb_diff) *nb_diff = diff;
    return diff == 0;
}

// ============================================================================
// PRIMITIVES D'INTERSECTION
// ============================================================================

// Premier t dans [lo, hi) avec a[t] >= key, par recherche exponentielle
// depuis lo (coût O(log(t - lo))).
static inline long long gallop_lb(const int* a, long long lo, long long hi,
                                  int key) {
    if (lo >= hi || a[lo] >= key) return lo;
    long long step = 1;
    while (lo + step < hi && a[lo + step] < key) step <<= 1;
    long long L = lo + step / 2 + 1, R = min(lo + step, hi - 1) + 1;
    // invariant : a[lo + step/2] < key ; borne sup exclusive R
    while (L < R) {
        long long mid = L + (R - L) / 2;
        if (a[mid] < key) L = mid + 1; else R = mid;
    }
    return L;
}

struct FingerStats {
    long long queries = 0;      // paires (i,j) traitées
    long long rewinds = 0;      // déclenchements du garde
    long long rewind_dist = 0;  // distance totale rembobinée
};

// Garde du Théorème 2 : position de départ sûre pour la requête courante.
// start = position sauvegardée ; [cb, ce) = colonne ; rmin = min du côté ligne.
static inline long long guarded_start(const int* cidx, long long cb,
                                      long long ce, long long start, int rmin,
                                      FingerStats& st) {
    st.queries++;
    if (start > cb && cidx[start - 1] >= rmin) {
        st.rewinds++;
        long long s = (long long)(lower_bound(cidx + cb, cidx + start, rmin)
                                  - cidx);
        st.rewind_dist += (start - s);
        return s;
    }
    return start;
}

// Intersection comptage (structure seule). gallop=true : avancée exponentielle.
// Retourne |a ∩ b| ; iy_end reçoit la position finale côté b (pour le finger).
static inline long long intersect_count(const int* a, long long ab, long long ae,
                                        const int* b, long long bb, long long be,
                                        long long iy_start, bool gallop,
                                        long long* iy_end = nullptr) {
    long long ix = ab, iy = iy_start, c = 0;
    while (ix < ae && iy < be) {
        int x = a[ix], y = b[iy];
        if (x == y) { c++; ix++; iy++; }
        else if (x < y) {
            ix = gallop ? gallop_lb(a, ix + 1, ae, y) : ix + 1;
        } else {
            iy = gallop ? gallop_lb(b, iy + 1, be, x) : iy + 1;
        }
    }
    if (iy_end) *iy_end = iy;
    return c;
}

// Intersection somme de produits (SpGEMM) : val alignées sur idx.
static inline double intersect_dot(const int* a, const double* av,
                                   long long ab, long long ae,
                                   const int* b, const double* bv,
                                   long long bb, long long be,
                                   long long iy_start, bool gallop,
                                   long long* iy_end = nullptr) {
    long long ix = ab, iy = iy_start;
    double s = 0;
    while (ix < ae && iy < be) {
        int x = a[ix], y = b[iy];
        if (x == y) { s += av[ix] * bv[iy]; ix++; iy++; }
        else if (x < y) {
            ix = gallop ? gallop_lb(a, ix + 1, ae, y) : ix + 1;
        } else {
            iy = gallop ? gallop_lb(b, iy + 1, be, x) : iy + 1;
        }
    }
    if (iy_end) *iy_end = iy;
    return s;
}

// ============================================================================
// SPGEMM COMPLET (section positionnement du papier)
// ============================================================================

// Référence : Gustavson avec SPA dense (accès positionnels purs).
static CSR spgemm_gustavson(const CSR& A, const CSR& B) {
    CSR C; C.n = A.n; C.m = B.m;
    C.indptr.assign((size_t)C.n + 1, 0);
    vector<double> spa((size_t)B.m, 0.0);
    vector<int> touched; touched.reserve(4096);
    vector<int> out_idx_row; vector<double> out_val_row;

    for (int i = 0; i < A.n; i++) {
        touched.clear();
        for (long long e = A.indptr[i]; e < A.indptr[(size_t)i + 1]; e++) {
            int k = A.idx[e]; double a = A.val[e];
            for (long long f = B.indptr[k]; f < B.indptr[(size_t)k + 1]; f++) {
                int j = B.idx[f];
                if (spa[j] == 0.0) touched.push_back(j);
                spa[j] += a * B.val[f];
            }
        }
        sort(touched.begin(), touched.end());
        for (int j : touched) {
            if (spa[j] != 0.0) { C.idx.push_back(j); C.val.push_back(spa[j]); }
            spa[j] = 0.0;
        }
        C.indptr[(size_t)i + 1] = C.nnz();
    }
    return C;
}

// Merge toutes-paires. mode: 0 = stateless, 1 = finger persistant NON GARDÉ
// (l'algorithme des chiffres TPDS -- conservé UNIQUEMENT pour la démonstration
// d'incorrection dans verify), 2 = finger persistant gardé (Théorème 2).
static CSR spgemm_merge(const CSR& A, const CSR& Bcsc, int mode, bool gallop,
                        FingerStats* stats = nullptr) {
    CSR C; C.n = A.n; C.m = Bcsc.n; // Bcsc = CSR(B^T) => Bcsc.n = colonnes de B
    C.indptr.assign((size_t)C.n + 1, 0);
    const int p = Bcsc.n;
    vector<long long> finger((size_t)p, -1); // -1 : réinitialisé par colonne? non:
    // fingers = position absolue dans Bcsc.idx ; init au début de chaque colonne
    for (int j = 0; j < p; j++) finger[j] = Bcsc.indptr[j];
    FingerStats st;

    for (int i = 0; i < A.n; i++) {
        long long rb = A.indptr[i], re = A.indptr[(size_t)i + 1];
        if (rb == re) { C.indptr[(size_t)i + 1] = C.nnz(); continue; }
        int rmin = A.idx[rb];
        for (int j = 0; j < p; j++) {
            long long cb = Bcsc.indptr[j], ce = Bcsc.indptr[(size_t)j + 1];
            if (cb == ce) continue;
            long long start;
            if (mode == 0) start = cb;
            else if (mode == 1) { st.queries++; start = finger[j]; }
            else start = guarded_start(Bcsc.idx.data(), cb, ce, finger[j],
                                       rmin, st);
            long long iy_end;
            double s = intersect_dot(A.idx.data(), A.val.data(), rb, re,
                                     Bcsc.idx.data(), Bcsc.val.data(), cb, ce,
                                     start, gallop, &iy_end);
            if (mode != 0) finger[j] = iy_end;
            if (s != 0.0) { C.idx.push_back(j); C.val.push_back(s); }
        }
        C.indptr[(size_t)i + 1] = C.nnz();
    }
    if (stats) *stats = st;
    return C;
}

// ============================================================================
// GRAPHES (comptage de triangles) : CSR symétrique, simple, sans boucle
// ============================================================================

static CSR make_graph(int n, vector<pair<int,int>>& edges) {
    vector<tuple<int,int,double>> t;
    t.reserve(edges.size() * 2);
    for (auto [u, v] : edges) {
        if (u == v) continue;
        t.push_back({u, v, 1.0});
        t.push_back({v, u, 1.0});
    }
    CSR G = from_triplets(n, n, t);
    for (double& x : G.val) x = 1.0; // doublons fusionnés -> re-binariser
    return G;
}

static CSR gen_er(int n, double avg_deg, int seed = 42) {
    mt19937_64 rng(seed);
    long long m = (long long)(n * avg_deg / 2);
    vector<pair<int,int>> e; e.reserve(m);
    uniform_int_distribution<int> d(0, n - 1);
    for (long long k = 0; k < m; k++) e.push_back({d(rng), d(rng)});
    return make_graph(n, e);
}

static CSR gen_ba(int n, int madd, int seed = 42) {
    mt19937_64 rng(seed);
    vector<pair<int,int>> e;
    vector<int> targets; // liste à répétition (attachement préférentiel)
    int m0 = madd + 1;
    for (int i = 0; i < m0; i++)
        for (int j = i + 1; j < m0; j++) {
            e.push_back({i, j}); targets.push_back(i); targets.push_back(j);
        }
    for (int v = m0; v < n; v++) {
        for (int k = 0; k < madd; k++) {
            int u = targets[uniform_int_distribution<size_t>(0, targets.size()-1)(rng)];
            e.push_back({u, v});
            targets.push_back(u); targets.push_back(v);
        }
    }
    return make_graph(n, e);
}

static CSR gen_band_graph(int n, int w) {
    vector<pair<int,int>> e;
    for (int i = 0; i < n; i++)
        for (int d = 1; d <= w && i + d < n; d++) e.push_back({i, i + d});
    return make_graph(n, e);
}

static CSR load_snap(const string& file) {
    ifstream f(file);
    if (!f.is_open()) throw runtime_error("Impossible d'ouvrir: " + file);
    vector<pair<int,int>> e;
    vector<long long> raw_u, raw_v;
    string line;
    unordered_map<long long,int> remap;
    int next_id = 0;
    auto id_of = [&](long long x) {
        auto it = remap.find(x);
        if (it != remap.end()) return it->second;
        remap[x] = next_id; return next_id++;
    };
    while (getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '%') continue;
        istringstream iss(line);
        long long u, v;
        if (!(iss >> u >> v)) continue;
        e.push_back({id_of(u), id_of(v)});
    }
    return make_graph(next_id, e);
}

static CSR load_mtx(const string& file) {
    ifstream f(file);
    if (!f.is_open()) throw runtime_error("Impossible d'ouvrir: " + file);
    string line;
    getline(f, line); // header %%MatrixMarket
    bool sym = line.find("symmetric") != string::npos;
    while (getline(f, line) && !line.empty() && line[0] == '%') {}
    istringstream iss(line);
    long long nr, nc, nnz;
    iss >> nr >> nc >> nnz;
    vector<pair<int,int>> e; e.reserve(nnz);
    for (long long k = 0; k < nnz; k++) {
        if (!getline(f, line)) break;
        istringstream ls(line);
        long long i, j; double v;
        ls >> i >> j; // valeur éventuelle ignorée (structure)
        (void)v;
        e.push_back({(int)(i - 1), (int)(j - 1)});
        if (!sym) {} // make_graph symétrise de toute façon (graphe non orienté)
    }
    return make_graph((int)max(nr, nc), e);
}

static CSR load_graph_arg(const string& arg) {
    if (arg.rfind("gen:", 0) == 0) {
        vector<string> parts;
        stringstream ss(arg);
        string tok;
        while (getline(ss, tok, ':')) parts.push_back(tok);
        if (parts.size() >= 4 && parts[1] == "er")
            return gen_er(stoi(parts[2]), stod(parts[3]));
        if (parts.size() >= 4 && parts[1] == "ba")
            return gen_ba(stoi(parts[2]), stoi(parts[3]));
        if (parts.size() >= 4 && parts[1] == "band")
            return gen_band_graph(stoi(parts[2]), stoi(parts[3]));
        throw runtime_error("Générateur inconnu: " + arg);
    }
    if (arg.size() > 4 && arg.substr(arg.size() - 4) == ".mtx")
        return load_mtx(arg);
    return load_snap(arg);
}

// ============================================================================
// ORDONNANCEMENTS (contribution C : le rewind rate dépend de l'ordre)
// ============================================================================

static vector<int> order_rcm(const CSR& G) {
    int n = G.n;
    vector<int> deg(n);
    for (int i = 0; i < n; i++)
        deg[i] = (int)(G.indptr[(size_t)i + 1] - G.indptr[i]);
    vector<char> vis(n, 0);
    vector<int> order; order.reserve(n);
    vector<int> nodes(n);
    iota(nodes.begin(), nodes.end(), 0);
    sort(nodes.begin(), nodes.end(),
         [&](int a, int b) { return deg[a] < deg[b]; });
    vector<int> nb;
    for (int s : nodes) {
        if (vis[s]) continue;
        queue<int> q; q.push(s); vis[s] = 1;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            order.push_back(u);
            nb.clear();
            for (long long e = G.indptr[u]; e < G.indptr[(size_t)u + 1]; e++)
                if (!vis[G.idx[e]]) nb.push_back(G.idx[e]);
            sort(nb.begin(), nb.end(),
                 [&](int a, int b) { return deg[a] < deg[b]; });
            for (int v : nb) if (!vis[v]) { vis[v] = 1; q.push(v); }
        }
    }
    reverse(order.begin(), order.end()); // Reverse Cuthill-McKee
    vector<int> old2new(n);
    for (int k = 0; k < n; k++) old2new[order[k]] = k;
    return old2new;
}

static vector<int> order_degree(const CSR& G) {
    int n = G.n;
    vector<int> nodes(n);
    iota(nodes.begin(), nodes.end(), 0);
    sort(nodes.begin(), nodes.end(), [&](int a, int b) {
        long long da = G.indptr[(size_t)a+1] - G.indptr[a];
        long long db = G.indptr[(size_t)b+1] - G.indptr[b];
        return make_pair(da, (long long)a) < make_pair(db, (long long)b);
    });
    vector<int> old2new(n);
    for (int k = 0; k < n; k++) old2new[nodes[k]] = k;
    return old2new;
}

static CSR permute_graph(const CSR& G, const vector<int>& old2new) {
    vector<pair<int,int>> e;
    e.reserve((size_t)(G.nnz() / 2));
    for (int i = 0; i < G.n; i++)
        for (long long f = G.indptr[i]; f < G.indptr[(size_t)i + 1]; f++) {
            int j = G.idx[f];
            if (i < j) e.push_back({old2new[i], old2new[j]});
        }
    return make_graph(G.n, e);
}

// ============================================================================
// COMPTAGE DE TRIANGLES (contribution phare)
// Sémantique : pour chaque arête (i,j) i<j, |N(i) ∩ N(j)| ; total = 3T.
// ============================================================================

// Gustavson masqué : ligne complète de A^2 puis filtrage par le masque.
// Coût = flops de A^2 -> montre l'inadéquation de Gustavson sous masque.
static long long tri_masked_gustavson(const CSR& G) {
    vector<long long> spa((size_t)G.n, 0);
    vector<int> touched;
    long long total = 0;
    for (int i = 0; i < G.n; i++) {
        touched.clear();
        for (long long e = G.indptr[i]; e < G.indptr[(size_t)i + 1]; e++) {
            int k = G.idx[e];
            for (long long f = G.indptr[k]; f < G.indptr[(size_t)k + 1]; f++) {
                int j = G.idx[f];
                if (spa[j] == 0) touched.push_back(j);
                spa[j]++;
            }
        }
        // application du masque : somme sur j ∈ N(i)
        for (long long e = G.indptr[i]; e < G.indptr[(size_t)i + 1]; e++)
            total += spa[G.idx[e]];
        for (int j : touched) spa[j] = 0;
    }
    return total / 6; // chaque triangle compté 2x par sommet, 3 sommets
}

// Edge-iterator : mode 0 = stateless, 2 = finger gardé ; gallop au choix.
// T threads : blocs de lignes contigus, fingers privés par thread
// (la théorie s'applique par bloc, i croissant dans chaque bloc).
static long long tri_edge(const CSR& G, int mode, bool gallop, int T,
                          FingerStats* stats = nullptr) {
    long long total = 0;
    FingerStats agg;
    omp_set_num_threads(T);
#pragma omp parallel reduction(+:total)
    {
        vector<long long> finger;
        if (mode == 2) {
            finger.resize((size_t)G.n);
            for (int j = 0; j < G.n; j++) finger[j] = G.indptr[j];
        }
        FingerStats st;
#pragma omp for schedule(runtime)
        for (int i = 0; i < G.n; i++) {
            long long rb = G.indptr[i], re = G.indptr[(size_t)i + 1];
            if (rb == re) continue;
            int rmin = G.idx[rb];
            for (long long e = rb; e < re; e++) {
                int j = G.idx[e];
                if (j <= i) continue; // arêtes i<j
                long long cb = G.indptr[j], ce = G.indptr[(size_t)j + 1];
                long long start = cb;
                if (mode == 2)
                    start = guarded_start(G.idx.data(), cb, ce, finger[j],
                                          rmin, st);
                long long iy_end;
                total += intersect_count(G.idx.data(), rb, re,
                                         G.idx.data(), cb, ce,
                                         start, gallop, &iy_end);
                if (mode == 2) finger[j] = iy_end;
            }
        }
#pragma omp critical
        {
            agg.queries += st.queries;
            agg.rewinds += st.rewinds;
            agg.rewind_dist += st.rewind_dist;
        }
    }
    if (stats) *stats = agg;
    return total / 3;
}

// Compact-forward : relabel par degré croissant, listes forward F(i) =
// {j ∈ N(i) : j > i} ; T = somme sur arêtes (i,j) i<j de |F(i) ∩ F(j)|.
// mode 0 = merge, gallop au choix, mode 2 = fingers persistants sur F(j).
struct ForwardGraph { CSR F; };

static ForwardGraph build_forward(const CSR& G_deg_ordered) {
    const CSR& G = G_deg_ordered;
    ForwardGraph FG;
    FG.F.n = G.n; FG.F.m = G.n;
    FG.F.indptr.assign((size_t)G.n + 1, 0);
    for (int i = 0; i < G.n; i++) {
        for (long long e = G.indptr[i]; e < G.indptr[(size_t)i + 1]; e++)
            if (G.idx[e] > i) {
                FG.F.idx.push_back(G.idx[e]);
                FG.F.val.push_back(1.0);
            }
        FG.F.indptr[(size_t)i + 1] = FG.F.nnz();
    }
    return FG;
}

static long long tri_cf(const CSR& G_deg_ordered, int mode, bool gallop, int T,
                        FingerStats* stats = nullptr) {
    const CSR& G = G_deg_ordered;
    ForwardGraph FG = build_forward(G);
    const CSR& F = FG.F;
    long long total = 0;
    FingerStats agg;
    omp_set_num_threads(T);
#pragma omp parallel reduction(+:total)
    {
        vector<long long> finger;
        if (mode == 2) {
            finger.resize((size_t)F.n);
            for (int j = 0; j < F.n; j++) finger[j] = F.indptr[j];
        }
        FingerStats st;
#pragma omp for schedule(runtime)
        for (int i = 0; i < F.n; i++) {
            long long rb = F.indptr[i], re = F.indptr[(size_t)i + 1];
            if (rb == re) continue;
            int rmin = F.idx[rb];
            for (long long e = rb; e < re; e++) {
                int j = F.idx[e]; // j > i par construction
                long long cb = F.indptr[j], ce = F.indptr[(size_t)j + 1];
                if (cb == ce) continue;
                long long start = cb;
                if (mode == 2)
                    start = guarded_start(F.idx.data(), cb, ce, finger[j],
                                          rmin, st);
                long long iy_end;
                total += intersect_count(F.idx.data(), rb, re,
                                         F.idx.data(), cb, ce,
                                         start, gallop, &iy_end);
                if (mode == 2) finger[j] = iy_end;
            }
        }
#pragma omp critical
        {
            agg.queries += st.queries;
            agg.rewinds += st.rewinds;
            agg.rewind_dist += st.rewind_dist;
        }
    }
    if (stats) *stats = agg;
    return total; // compté exactement une fois par triangle
}

// ============================================================================
// BENCH
// ============================================================================

// Dispersion de la derniere serie (min, mediane, max) -- pour la section
// methodologie exigee par les relecteurs : imprimee sur stdout apres chaque
// benchmark, sans changer le schema des CSV.
static double g_last_min = 0, g_last_max = 0;
static double bench_ms(const function<void()>& f, int runs) {
    vector<double> t;
    f(); // echauffement (non mesure)
    for (int r = 0; r < runs; r++) {
        auto t0 = chrono::high_resolution_clock::now();
        f();
        auto t1 = chrono::high_resolution_clock::now();
        t.push_back(chrono::duration<double, milli>(t1 - t0).count());
    }
    sort(t.begin(), t.end());
    g_last_min = t.front(); g_last_max = t.back();
    double med = t[t.size() / 2];
    cout << "    [spread] min=" << g_last_min << " med=" << med
         << " max=" << g_last_max
         << " cv~" << (med > 0 ? (g_last_max - g_last_min) / med : 0) << endl;
    return med;
}

// ============================================================================
// GÉNÉRATEURS SPGEMM COMPLET (beta / bande) -- version CSR directe
// ============================================================================

static CSR gen_beta_matrix(int n, double density, double beta, int seed) {
    mt19937 rng(seed);
    int per_col = max(1, (int)(n * density));
    vector<tuple<int,int,double>> t;
    vector<int> prev;
    uniform_int_distribution<int> d(0, n - 1);
    for (int k = 0; k < per_col; k++) prev.push_back(d(rng));
    sort(prev.begin(), prev.end());
    prev.erase(unique(prev.begin(), prev.end()), prev.end());
    for (int r : prev) t.push_back({r, 0, 1.0});
    for (int j = 1; j < n; j++) {
        vector<int> curr;
        int keep = (int)(beta * prev.size());
        shuffle(prev.begin(), prev.end(), rng);
        for (int k = 0; k < keep && k < (int)prev.size(); k++)
            curr.push_back(prev[k]);
        vector<char> in(n, 0);
        for (int x : curr) in[x] = 1;
        while ((int)curr.size() < (int)prev.size() && (int)curr.size() < n) {
            int r = d(rng);
            if (!in[r]) { in[r] = 1; curr.push_back(r); }
        }
        sort(curr.begin(), curr.end());
        curr.erase(unique(curr.begin(), curr.end()), curr.end());
        for (int r : curr) t.push_back({r, j, 1.0});
        prev = curr;
    }
    return from_triplets(n, n, t);
}

static CSR gen_band_matrix(int n, int w) {
    vector<tuple<int,int,double>> t;
    for (int i = 0; i < n; i++)
        for (int k = i; k < min(n, i + w); k++) t.push_back({i, k, 1.0});
    return from_triplets(n, n, t);
}

// ============================================================================
// EXPÉRIENCES
// ============================================================================

static void experiment_verify() {
    cout << "=== VERIFY ===\n(1 = identique a la reference, 0 = DIVERGENT)\n\n";

    // ---- SpGEMM complet ----
    cout << "-- SpGEMM complet (reference: Gustavson SPA) --" << endl;
    struct Cfg { const char* name; CSR A, B; };
    vector<Cfg> cfgs;
    cfgs.push_back({"random n=512 d=0.02 b=0.6",
                    gen_beta_matrix(512, 0.02, 0.6, 42),
                    gen_beta_matrix(512, 0.02, 0.6, 1000)});
    cfgs.push_back({"band n=1024 w=16",
                    gen_band_matrix(1024, 16), gen_band_matrix(1024, 16)});
    cfgs.push_back({"random n=2048 d=0.01 b=0.6",
                    gen_beta_matrix(2048, 0.01, 0.6, 42),
                    gen_beta_matrix(2048, 0.01, 0.6, 1000)});

    for (auto& c : cfgs) {
        CSR Bcsc = transpose(c.B);
        CSR ref = spgemm_gustavson(c.A, c.B);
        long long d0, d1, d2, d3;
        bool ok0 = same_csr(spgemm_merge(c.A, Bcsc, 0, false), ref, &d0);
        bool ok1 = same_csr(spgemm_merge(c.A, Bcsc, 1, false), ref, &d1);
        bool ok2 = same_csr(spgemm_merge(c.A, Bcsc, 2, false), ref, &d2);
        bool ok3 = same_csr(spgemm_merge(c.A, Bcsc, 2, true),  ref, &d3);
        // survie du non garde : nnz produits / nnz de la reference
        CSR wrong = spgemm_merge(c.A, Bcsc, 1, false);
        double surv = ref.nnz() ? (double)wrong.nnz() / ref.nnz() : 0.0;
        cout << left << setw(28) << c.name
             << " stateless=" << ok0
             << " NON_GARDE=" << ok1 << "(diff=" << d1
             << ", survie=" << surv * 100 << "%)"
             << " garde=" << ok2
             << " garde+gallop=" << ok3 << endl;
    }

    // ---- Triangles ----
    cout << "\n-- Triangles (cross-check: Gustavson masque vs edge-merge) --" << endl;
    vector<pair<string, CSR>> graphs;
    graphs.push_back({"er n=2000 deg=20", gen_er(2000, 20)});
    graphs.push_back({"ba n=2000 m=8", gen_ba(2000, 8)});
    graphs.push_back({"band n=2000 w=12", gen_band_graph(2000, 12)});

    for (auto& [name, G] : graphs) {
        long long t_ref = tri_masked_gustavson(G);
        long long t_m  = tri_edge(G, 0, false, 1);
        long long t_g  = tri_edge(G, 0, true, 1);
        long long t_f  = tri_edge(G, 2, true, 1);
        long long t_fp = tri_edge(G, 2, true, 4);
        CSR Gd = permute_graph(G, order_degree(G));
        long long t_cf  = tri_cf(Gd, 0, true, 1);
        long long t_cff = tri_cf(Gd, 2, true, 1);
        bool all_ok = (t_m == t_ref && t_g == t_ref && t_f == t_ref &&
                       t_fp == t_ref && t_cf == t_ref && t_cff == t_ref);
        cout << left << setw(22) << name << " T=" << setw(10) << t_ref
             << " merge=" << (t_m == t_ref) << " gallop=" << (t_g == t_ref)
             << " finger=" << (t_f == t_ref) << " finger_par=" << (t_fp == t_ref)
             << " cf=" << (t_cf == t_ref) << " cf_finger=" << (t_cff == t_ref)
             << (all_ok ? "  [OK]" : "  [ECHEC]") << endl;
    }
}

static void experiment_fullspgemm(const string& out, int n, int trials,
                                  int runs) {
    cout << "=== FULLSPGEMM (positionnement) ===" << endl;
    ofstream csv(out);
    csv << "kind,param,n,nnz,rho_A,rho_B,"
        << "t_gustavson_ms,t_stateless_ms,t_stateless_gallop_ms,"
        << "t_guarded_ms,t_guarded_gallop_ms,"
        << "ratio_stateless_over_guarded,ratio_stateless_over_guarded_gallop,"
        << "bound_1_plus_rhoB_over_rhoA,rewind_rate,"
        << "correct_guarded_gallop" << endl;

    auto run_case = [&](const string& kind, double param, CSR& A, CSR& B) {
        CSR Bcsc = transpose(B);
        double rho_A = (double)A.nnz() / A.n;
        double rho_B = (double)B.nnz() / B.m;

        CSR ref = spgemm_gustavson(A, B);
        long long dd;
        bool ok = same_csr(spgemm_merge(A, Bcsc, 2, true), ref, &dd);

        double t_gus = bench_ms([&]{ auto C = spgemm_gustavson(A, B); g_sink += C.nnz(); }, runs);
        double t_sl  = bench_ms([&]{ auto C = spgemm_merge(A, Bcsc, 0, false); g_sink += C.nnz(); }, runs);
        double t_slg = bench_ms([&]{ auto C = spgemm_merge(A, Bcsc, 0, true); g_sink += C.nnz(); }, runs);
        double t_gd  = bench_ms([&]{ auto C = spgemm_merge(A, Bcsc, 2, false); g_sink += C.nnz(); }, runs);
        double t_gdg = bench_ms([&]{ auto C = spgemm_merge(A, Bcsc, 2, true); g_sink += C.nnz(); }, runs);

        FingerStats st;
        { auto C = spgemm_merge(A, Bcsc, 2, true, &st); (void)C; }
        double rr = st.queries ? (double)st.rewinds / st.queries : 0;

        csv << kind << "," << param << "," << A.n << "," << A.nnz() << ","
            << rho_A << "," << rho_B << ","
            << t_gus << "," << t_sl << "," << t_slg << ","
            << t_gd << "," << t_gdg << ","
            << (t_sl / t_gd) << "," << (t_sl / t_gdg) << ","
            << (1.0 + rho_B / max(1.0, rho_A)) << "," << rr << ","
            << ok << endl;
        cout << kind << " " << param
             << " | gus=" << t_gus << " sl=" << t_sl << " sl+g=" << t_slg
             << " gd=" << t_gd << " gd+g=" << t_gdg
             << " | rewind=" << rr << " correct=" << ok << endl;
    };

    for (double beta : {0.0, 0.3, 0.6, 0.9, 0.95}) {
        for (int tr = 0; tr < trials; tr++) {
            CSR A = gen_beta_matrix(n, 0.01, beta, 42 + tr);
            CSR B = gen_beta_matrix(n, 0.01, beta, 1000 + tr);
            run_case("beta", beta, A, B);
        }
    }
    // 128-512 : question du croisement avec Gustavson (section 7.5 de
    // l'article) ; sur Raspberry Pi, s'arreter a 64 est acceptable.
    for (int w : {16, 32, 64, 128, 256, 512}) {
        CSR A = gen_band_matrix(n, w);
        CSR B = gen_band_matrix(n, w);
        run_case("band", w, A, B);
    }
    csv.close();
    cout << "Sauvegarde: " << out << endl;
}

static void experiment_triangles(const string& out,
                                 const vector<string>& graph_args,
                                 int runs, int T) {
    cout << "=== TRIANGLES ===" << endl;
    ofstream csv(out);
    csv << "graph,ordering,n,nnz,avg_deg,triangles,"
        << "t_masked_gustavson_ms,t_edge_merge_ms,t_edge_gallop_ms,"
        << "t_edge_finger_ms,t_cf_gallop_ms,t_cf_finger_ms,"
        << "rewind_rate_edge_finger,rewind_rate_cf_finger,"
        << "t_edge_finger_par_ms,t_cf_finger_par_ms,threads,"
        << "counts_all_equal" << endl;

    for (const string& arg : graph_args) {
        cout << "\nGraphe: " << arg << endl;
        CSR G0 = load_graph_arg(arg);
        cout << "  n=" << G0.n << " nnz=" << G0.nnz()
             << " deg_moy=" << (double)G0.nnz() / G0.n << endl;

        struct Ord { const char* name; CSR G; };
        vector<Ord> ords;
        ords.push_back({"natural", G0});
        ords.push_back({"rcm", permute_graph(G0, order_rcm(G0))});
        ords.push_back({"degree", permute_graph(G0, order_degree(G0))});

        for (auto& o : ords) {
            const CSR& G = o.G;
            CSR Gd = permute_graph(G, order_degree(G)); // pour compact-forward

            long long c_ref = tri_masked_gustavson(G);
            long long c_m  = tri_edge(G, 0, false, 1);
            long long c_g  = tri_edge(G, 0, true, 1);
            FingerStats st_e, st_cf;
            long long c_f  = tri_edge(G, 2, true, 1, &st_e);
            long long c_cg = tri_cf(Gd, 0, true, 1);
            long long c_cf = tri_cf(Gd, 2, true, 1, &st_cf);
            bool eq = (c_m == c_ref && c_g == c_ref && c_f == c_ref &&
                       c_cg == c_ref && c_cf == c_ref);

            double t_mg = bench_ms([&]{ g_sink += tri_masked_gustavson(G); }, runs);
            double t_m  = bench_ms([&]{ g_sink += tri_edge(G, 0, false, 1); }, runs);
            double t_g  = bench_ms([&]{ g_sink += tri_edge(G, 0, true, 1); }, runs);
            double t_f  = bench_ms([&]{ g_sink += tri_edge(G, 2, true, 1); }, runs);
            double t_cg = bench_ms([&]{ g_sink += tri_cf(Gd, 0, true, 1); }, runs);
            double t_cf = bench_ms([&]{ g_sink += tri_cf(Gd, 2, true, 1); }, runs);
            double t_fp = (T > 1)
                ? bench_ms([&]{ g_sink += tri_edge(G, 2, true, T); }, runs) : -1;
            double t_cp = (T > 1)
                ? bench_ms([&]{ g_sink += tri_cf(Gd, 2, true, T); }, runs) : -1;

            double rr_e = st_e.queries ? (double)st_e.rewinds / st_e.queries : 0;
            double rr_c = st_cf.queries ? (double)st_cf.rewinds / st_cf.queries : 0;

            csv << arg << "," << o.name << "," << G.n << "," << G.nnz() << ","
                << (double)G.nnz() / G.n << "," << c_ref << ","
                << t_mg << "," << t_m << "," << t_g << "," << t_f << ","
                << t_cg << "," << t_cf << ","
                << rr_e << "," << rr_c << ","
                << t_fp << "," << t_cp << "," << T << ","
                << eq << endl;

            cout << "  [" << o.name << "] T=" << c_ref
                 << " | mgus=" << t_mg << " merge=" << t_m << " gallop=" << t_g
                 << " finger=" << t_f << " (rr=" << rr_e << ")"
                 << " | cf=" << t_cg << " cf_finger=" << t_cf
                 << " (rr=" << rr_c << ")"
                 << " | counts_ok=" << eq << endl;
        }
    }
    csv.close();
    cout << "\nSauvegarde: " << out << endl;
}

// Sélecteur adaptatif : échantillonne P% des lignes avec le noyau finger,
// mesure le rewind rate, choisit finger si rate < seuil sinon gallop pur.
static void experiment_adaptive(const string& out,
                                const vector<string>& graph_args,
                                double sample_pct, double threshold,
                                int runs) {
    cout << "=== ADAPTIVE (rewind rate comme predicteur) ===" << endl;
    ofstream csv(out);
    csv << "graph,sample_pct,sampled_rewind_rate,choice,threshold,"
        << "t_sampling_ms,t_chosen_ms,t_total_adaptive_ms,"
        << "t_gallop_ms,t_finger_ms,t_best_oracle_ms,overhead_vs_oracle"
        << endl;

    for (const string& arg : graph_args) {
        CSR G = load_graph_arg(arg);
        int n_sample = max(1, (int)(G.n * sample_pct / 100.0));

        // échantillonnage : finger sur les n_sample premières lignes
        FingerStats st;
        auto t0 = chrono::high_resolution_clock::now();
        {
            vector<long long> finger((size_t)G.n);
            for (int j = 0; j < G.n; j++) finger[j] = G.indptr[j];
            for (int i = 0; i < n_sample; i++) {
                long long rb = G.indptr[i], re = G.indptr[(size_t)i + 1];
                if (rb == re) continue;
                int rmin = G.idx[rb];
                for (long long e = rb; e < re; e++) {
                    int j = G.idx[e];
                    if (j <= i) continue;
                    long long cb = G.indptr[j], ce = G.indptr[(size_t)j + 1];
                    long long start = guarded_start(G.idx.data(), cb, ce,
                                                    finger[j], rmin, st);
                    long long iy_end;
                    intersect_count(G.idx.data(), rb, re, G.idx.data(), cb, ce,
                                    start, true, &iy_end);
                    finger[j] = iy_end;
                }
            }
        }
        auto t1 = chrono::high_resolution_clock::now();
        double t_samp = chrono::duration<double, milli>(t1 - t0).count();
        double rr = st.queries ? (double)st.rewinds / st.queries : 0;
        bool choose_finger = rr < threshold;

        double t_gal = bench_ms([&]{ g_sink += tri_edge(G, 0, true, 1); }, runs);
        double t_fin = bench_ms([&]{ g_sink += tri_edge(G, 2, true, 1); }, runs);
        double t_chosen = choose_finger ? t_fin : t_gal;
        double t_oracle = min(t_gal, t_fin);

        csv << arg << "," << sample_pct << "," << rr << ","
            << (choose_finger ? "finger" : "gallop") << "," << threshold << ","
            << t_samp << "," << t_chosen << "," << (t_samp + t_chosen) << ","
            << t_gal << "," << t_fin << "," << t_oracle << ","
            << ((t_samp + t_chosen) / t_oracle) << endl;

        cout << arg << " | rr_echantillon=" << rr
             << " -> " << (choose_finger ? "finger" : "gallop")
             << " | gallop=" << t_gal << " finger=" << t_fin
             << " adaptatif_total=" << (t_samp + t_chosen) << endl;
    }
    csv.close();
    cout << "Sauvegarde: " << out << endl;
}

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage:\n"
             << "  " << argv[0] << " verify\n"
             << "  " << argv[0] << " fullspgemm <out.csv> [n=2048] [trials=3] [runs=3]\n"
             << "  " << argv[0] << " triangles <out.csv> [--runs R] [--threads T] <graphe...>\n"
             << "  " << argv[0] << " adaptive  <out.csv> [--sample-pct P] [--threshold S] <graphe...>\n"
             << "graphe: fichier.txt | fichier.mtx | gen:er:n:deg | gen:ba:n:m | gen:band:n:w\n";
        return 1;
    }
    string cmd = argv[1];
    try {
        if (cmd == "verify") {
            experiment_verify();
        } else if (cmd == "fullspgemm") {
            if (argc < 3) throw runtime_error("fullspgemm: <out.csv> requis");
            int n = argc >= 4 ? stoi(argv[3]) : 2048;
            int trials = argc >= 5 ? stoi(argv[4]) : 3;
            int runs = argc >= 6 ? stoi(argv[5]) : 3;
            experiment_fullspgemm(argv[2], n, trials, runs);
        } else if (cmd == "triangles") {
            if (argc < 4) throw runtime_error("triangles: <out.csv> <graphe...> requis");
            int runs = 3, T = 1;
            vector<string> graphs;
            for (int i = 3; i < argc; i++) {
                string a = argv[i];
                if (a == "--runs" && i + 1 < argc) runs = stoi(argv[++i]);
                else if (a == "--threads" && i + 1 < argc) T = stoi(argv[++i]);
                else graphs.push_back(a);
            }
            experiment_triangles(argv[2], graphs, runs, T);
        } else if (cmd == "adaptive") {
            if (argc < 4) throw runtime_error("adaptive: <out.csv> <graphe...> requis");
            double pct = 2.0, thr = 0.05;
            int runs = 3;
            vector<string> graphs;
            for (int i = 3; i < argc; i++) {
                string a = argv[i];
                if (a == "--sample-pct" && i + 1 < argc) pct = stod(argv[++i]);
                else if (a == "--threshold" && i + 1 < argc) thr = stod(argv[++i]);
                else if (a == "--runs" && i + 1 < argc) runs = stoi(argv[++i]);
                else graphs.push_back(a);
            }
            experiment_adaptive(argv[2], graphs, pct, thr, runs);
        } else {
            throw runtime_error("commande inconnue: " + cmd);
        }
    } catch (const exception& e) {
        cerr << "ERREUR: " << e.what() << endl;
        return 1;
    }
    return 0;
}
