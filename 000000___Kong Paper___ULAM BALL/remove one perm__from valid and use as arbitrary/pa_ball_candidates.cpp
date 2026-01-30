// pa_ball_candidates.cpp
// 1) Read PA file (one permutation per line).
// 2) Pick random center perm π from PA; remove it from PA file (overwrite).
// 3) Compute Ulam ball B(π, r=d-1) by translating B(id,r) built via RSK+InverseRSK.
// 4) For each τ in B(π,r), test if dU(τ, sigma) >= d for all sigma in PA.
// 5) Output candidates and counts.
//
// Build (WSL):
//   g++ -O2 -std=c++17 pa_ball_candidates.cpp -o pa_ball_candidates
//
// Run example (n=11, d=3):
//   ./pa_ball_candidates 11 3 11_3_pa.txt 11_3_ball.txt 11_3_candidates.txt
//
// Optional seed (repeatable randomness):
//   ./pa_ball_candidates 11 3 11_3_pa.txt 11_3_ball.txt 11_3_candidates.txt --seed 12345

#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using boost::multiprecision::cpp_int;

using Perm = vector<int>;
using Tableau = vector<vector<int>>;

// -------------------- file I/O --------------------
static bool read_perm_line(const string& line, int n, Perm& out){
    out.clear();
    out.reserve(n);
    stringstream ss(line);
    int x;
    while(ss >> x) out.push_back(x);
    return (int)out.size() == n;
}

static vector<Perm> read_perm_file(const string& filename, int n){
    ifstream fin(filename);
    if(!fin) throw runtime_error("Cannot open file: " + filename);

    vector<Perm> perms;
    string line;
    while(getline(fin, line)){
        if(line.empty()) continue;
        Perm p;
        if(!read_perm_line(line, n, p))
            throw runtime_error("Bad line (wrong length) in file: " + filename + "\nLine: " + line);
        perms.push_back(std::move(p));
    }
    return perms;
}

static void write_perm_file(const string& filename, const vector<Perm>& perms){
    ofstream fout(filename);
    if(!fout) throw runtime_error("Cannot write file: " + filename);
    for(const auto& p : perms){
        for(size_t i=0;i<p.size();i++){
            fout << p[i] << (i+1<p.size() ? ' ' : '\n');
        }
    }
}

static bool is_valid_perm_1_to_n(const Perm& p){
    int n = (int)p.size();
    vector<int> seen(n+1,0);
    for(int x: p){
        if(x<1 || x>n) return false;
        if(seen[x]) return false;
        seen[x]=1;
    }
    return true;
}

// -------------------- LIS with early stop --------------------
static int lis_len_early(const vector<int>& seq, int stop_at){
    // Returns LIS length, but stops and returns >=stop_at as soon as it reaches stop_at.
    vector<int> tail;
    tail.reserve(seq.size());
    for(int x: seq){
        auto it = lower_bound(tail.begin(), tail.end(), x);
        if(it==tail.end()) tail.push_back(x);
        else *it = x;
        if((int)tail.size() >= stop_at) return (int)tail.size();
    }
    return (int)tail.size();
}

// Ulam distance between permutations a and b using LIS of position sequence of a in b
static int ulam_distance(const Perm& a, const Perm& b, int early_fail_d){
    // early_fail_d = d: we can early stop if we already know distance < d.
    // distance < d  <=>  LCS > n-d  <=> LCS >= n-d+1
    int n = (int)a.size();
    int needLCS = n - early_fail_d + 1; // if LCS reaches this, then dist <= d-1 => FAIL

    // pos in b
    vector<int> pos(n+1);
    for(int i=0;i<n;i++) pos[b[i]] = i+1; // 1..n

    vector<int> seq(n);
    for(int i=0;i<n;i++) seq[i] = pos[a[i]];

    int L = lis_len_early(seq, needLCS);
    int dist = n - L;
    return dist;
}

// -------------------- Partitions for shapes λ with λ1 >= minFirst --------------------
static void gen_partitions_rec(int remaining, int maxPart, vector<int>& cur,
                               int minFirst, vector<vector<int>>& out){
    if(remaining==0){
        if(!cur.empty() && cur[0] >= minFirst) out.push_back(cur);
        return;
    }
    for(int x=min(maxPart, remaining); x>=1; --x){
        cur.push_back(x);
        gen_partitions_rec(remaining-x, x, cur, minFirst, out);
        cur.pop_back();
    }
}

static vector<vector<int>> gen_partitions(int n, int minFirst){
    vector<vector<int>> out;
    vector<int> cur;
    gen_partitions_rec(n, n, cur, minFirst, out);
    sort(out.begin(), out.end(), [](const auto& a, const auto& b){
        if(a[0]!=b[0]) return a[0]>b[0];
        return a>b;
    });
    return out;
}

// -------------------- SYT generation by removing corners --------------------
static vector<int> corners_of_shape(const vector<int>& lam){
    vector<int> idx;
    for(int i=0;i<(int)lam.size();i++){
        if(i==(int)lam.size()-1 || lam[i] > lam[i+1]) idx.push_back(i);
    }
    return idx;
}

static vector<int> shape_remove_corner(const vector<int>& lam, int cornerRow){
    vector<int> lam2 = lam;
    lam2[cornerRow]--;
    if(lam2[cornerRow]==0) lam2.erase(lam2.begin()+cornerRow);
    return lam2;
}

static vector<Tableau> gen_SYT_from_shape_rec(const vector<int>& lam, int k){
    if(k==0) return vector<Tableau>{Tableau{}};

    vector<Tableau> result;
    auto corners = corners_of_shape(lam);

    for(int cr : corners){
        vector<int> lam2 = shape_remove_corner(lam, cr);
        auto sub = gen_SYT_from_shape_rec(lam2, k-1);

        bool rowRemoved = (cr < (int)lam.size() && lam[cr]==1);

        for(auto T : sub){
            if(rowRemoved) T.insert(T.begin()+cr, vector<int>{k});
            else          T[cr].push_back(k);
            result.push_back(std::move(T));
        }
    }
    return result;
}

static vector<Tableau> AllSYT(const vector<int>& lam){
    int n=0; for(int x: lam) n+=x;
    return gen_SYT_from_shape_rec(lam, n);
}

// -------------------- Inverse RSK (reverse row insertion) --------------------
static Perm InverseRSK(Tableau P, Tableau Q){
    int n=0; for(auto &row: Q) n += (int)row.size();
    vector<int> w(n+1, 0);

    for(int k=n; k>=1; --k){
        int r=-1, c=-1;
        for(int i=0;i<(int)Q.size();i++){
            for(int j=0;j<(int)Q[i].size();j++){
                if(Q[i][j]==k){ r=i; c=j; break; }
            }
            if(r!=-1) break;
        }
        if(r==-1) throw runtime_error("InverseRSK: cannot find k in Q");

        Q[r].erase(Q[r].begin()+c);
        if(Q[r].empty()) Q.erase(Q.begin()+r);

        int x = P[r][c];
        P[r].erase(P[r].begin()+c);
        if(P[r].empty()) P.erase(P.begin()+r);

        for(int i=r-1; i>=0; --i){
            auto &row = P[i];
            auto it = lower_bound(row.begin(), row.end(), x);
            if(it==row.begin()) throw runtime_error("InverseRSK: reverse bump failed");
            --it;
            int y = *it;
            *it = x;
            x = y;
        }
        w[k] = x;
    }

    Perm perm(n);
    for(int i=1;i<=n;i++) perm[i-1]=w[i];
    return perm;
}

// σ = π ∘ τ in one-line notation: σ_i = π_{τ_i}
static Perm compose_left(const Perm& pi, const Perm& tau){
    int n=(int)pi.size();
    Perm sigma(n);
    for(int i=0;i<n;i++) sigma[i] = pi[tau[i]-1];
    return sigma;
}

// Build B(id,r) exactly using RSK layers with λ1 >= n-r
static vector<Perm> build_ball_id(int n, int d){
    int r = d-1;
    int minFirst = n - r;

    auto shapes = gen_partitions(n, minFirst);
    vector<Perm> ball_id;

    for(const auto& lam : shapes){
        auto T = AllSYT(lam);
        ball_id.reserve(ball_id.size() + T.size()*T.size());
        for(const auto& P : T){
            for(const auto& Q : T){
                ball_id.push_back(InverseRSK(P,Q));
            }
        }
    }
    return ball_id;
}

// Translate ball: B(pi,r) = { pi ∘ tau : tau ∈ B(id,r) }
static vector<Perm> translate_ball(const Perm& pi, const vector<Perm>& ball_id){
    vector<Perm> ball_pi;
    ball_pi.reserve(ball_id.size());
    for(const auto& tau : ball_id){
        ball_pi.push_back(compose_left(pi, tau));
    }
    return ball_pi;
}

// -------------------- main workflow --------------------
int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if(argc < 6){
        cerr << "Usage:\n"
             << "  " << argv[0] << " n d pa_file ball_file candidate_file [--seed SEED]\n\n"
             << "Example:\n"
             << "  " << argv[0] << " 11 3 11_3_pa.txt 11_3_ball.txt 11_3_candidates.txt\n";
        return 1;
    }

    int n = stoi(argv[1]);
    int d = stoi(argv[2]);
    string pa_file = argv[3];
    string ball_file = argv[4];
    string cand_file = argv[5];

    bool use_seed=false;
    uint64_t seed=0;
    for(int i=6;i+1<argc;i++){
        if(string(argv[i])=="--seed"){
            use_seed=true;
            seed = stoull(argv[i+1]);
        }
    }

    // 1) Read PA
    auto pa = read_perm_file(pa_file, n);
    if(pa.empty()){
        cerr << "ERROR: PA file is empty.\n";
        return 1;
    }
    for(const auto& p : pa){
        if(!is_valid_perm_1_to_n(p)){
            cerr << "ERROR: PA file contains an invalid permutation.\n";
            return 1;
        }
    }

    // 2) Pick random center and remove it
    mt19937_64 rng;
    if(use_seed) rng.seed(seed);
    else rng.seed((uint64_t)chrono::high_resolution_clock::now().time_since_epoch().count());

    uniform_int_distribution<size_t> dist_idx(0, pa.size()-1);
    size_t idx = dist_idx(rng);

    Perm center = pa[idx];
    pa.erase(pa.begin()+idx);

    // overwrite PA file with one less permutation
    write_perm_file(pa_file, pa);

    // also save removed center (for record)
    {
        ofstream out("removed_center.txt");
        for(int i=0;i<n;i++) out << center[i] << (i+1<n?' ':'\n');
    }

    cout << "Selected random center (removed from PA): saved to removed_center.txt\n";
    cout << "Updated PA file written back: " << pa_file << "\n";
    cout << "Updated PA size = " << pa.size() << "\n";

    // 3) Build B(id,r) then translate to B(center,r)
    cout << "Building B(id, r=d-1) using RSK+InverseRSK...\n";
    auto ball_id = build_ball_id(n, d);
    cout << "|B(id, d-1)| = " << ball_id.size() << "\n";

    cout << "Translating to ball centered at chosen permutation...\n";
    auto ball_center = translate_ball(center, ball_id);

    // write ball file
    write_perm_file(ball_file, ball_center);
    cout << "Ball written to: " << ball_file << " (size " << ball_center.size() << ")\n";

    // 4) Check each ball permutation vs all PA permutations
    cout << "Checking candidates: require dU(tau, sigma) >= d for all sigma in PA...\n";

    vector<Perm> candidates;
    candidates.reserve(ball_center.size());

    for(const auto& tau : ball_center){
        bool ok = true;
        for(const auto& sigma : pa){
            int du = ulam_distance(tau, sigma, d);
            if(du < d){
                ok = false;
                break;
            }
        }
        if(ok) candidates.push_back(tau);
    }

    // write candidates
    write_perm_file(cand_file, candidates);

    cout << "Candidates written to: " << cand_file << "\n";
    cout << "Total candidates found = " << candidates.size() << "\n";
    return 0;
}
