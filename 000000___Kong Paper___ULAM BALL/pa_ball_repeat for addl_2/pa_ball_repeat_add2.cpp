// pa_ball_repeat_add2.cpp
// Goal: k=1 augmentation for an (n,d)-PA under Ulam metric.
// Repeat trials:
//  1) pick random center pi from PA, remove it
//  2) compute ball B(pi, r=d-1) via translation from B(id,r)
//  3) candidates = elements of B(pi,r) that are distance >= d from all perms in PA\{pi}
//  4) if there exist two candidates c1,c2 with dU(c1,c2) >= d, then
//        PA <- (PA \ {pi}) U {c1,c2}
//     and write updated PA.
//
// Notes:
// - We compute B(id,r) once using RSK+InverseRSK construction (Kong style).
// - Ulam distance between permutations is computed as n - LIS(pos_in_sigma of pi's symbols).

#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using boost::multiprecision::cpp_int;

using Perm = vector<int>;
using Tableau = vector<vector<int>>;

// -------------------- factorial (big int) --------------------
static cpp_int factorial_int(int n){
    cpp_int f = 1;
    for(int i=2;i<=n;i++) f *= i;
    return f;
}

// -------------------- Hook-length formula f^lambda --------------------
static cpp_int f_lambda(const vector<int>& lam){
    int n=0; for(int x: lam) n += x;
    cpp_int num = factorial_int(n);
    cpp_int den = 1;

    for(int i=0;i<(int)lam.size();i++){
        for(int j=0;j<lam[i];j++){
            int right = lam[i] - (j+1);
            int below = 0;
            for(int k=i+1;k<(int)lam.size();k++){
                if(lam[k] > j) below++;
            }
            int hook = right + below + 1;
            den *= hook;
        }
    }
    return num/den;
}

// -------------------- partitions of n with first part >= minFirst --------------------
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
            if(rowRemoved){
                T.insert(T.begin()+cr, vector<int>{k});
            }else{
                T[cr].push_back(k);
            }
            result.push_back(std::move(T));
        }
    }
    return result;
}

static vector<Tableau> gen_SYT_from_shape(const vector<int>& lam){
    int n=0; for(int x: lam) n+=x;
    return gen_SYT_from_shape_rec(lam, n);
}

// -------------------- Inverse RSK (reverse row insertion) --------------------
static vector<int> inverse_rsk(Tableau P, Tableau Q){
    int n=0;
    for(auto &row: Q) n += (int)row.size();

    vector<int> w(n+1, 0); // 1-indexed

    for(int k=n; k>=1; --k){
        int r=-1, c=-1;
        for(int i=0;i<(int)Q.size();i++){
            for(int j=0;j<(int)Q[i].size();j++){
                if(Q[i][j]==k){ r=i; c=j; break; }
            }
            if(r!=-1) break;
        }
        if(r==-1) throw runtime_error("inverse_rsk: cannot find k in Q");

        Q[r].erase(Q[r].begin()+c);
        if(Q[r].empty()) Q.erase(Q.begin()+r);

        int x = P[r][c];
        P[r].erase(P[r].begin()+c);
        if(P[r].empty()) P.erase(P.begin()+r);

        // reverse bump upwards
        for(int i=r-1; i>=0; --i){
            auto &row = P[i];
            auto it = lower_bound(row.begin(), row.end(), x); // first >= x
            if(it==row.begin()){
                throw runtime_error("inverse_rsk: reverse bump failed");
            }
            --it; // rightmost < x
            int y = *it;
            *it = x;
            x = y;
        }
        w[k] = x;
    }

    vector<int> perm(n);
    for(int i=1;i<=n;i++) perm[i-1]=w[i];
    return perm;
}

// -------------------- compose (left translate): pi ∘ tau --------------------
static Perm compose_left(const Perm& pi, const Perm& tau){
    int n = (int)pi.size();
    Perm sigma(n);
    for(int i=0;i<n;i++){
        sigma[i] = pi[tau[i]-1];
    }
    return sigma;
}

// -------------------- LIS length --------------------
static int lis_len_positions(const vector<int>& seq){
    vector<int> tail;
    tail.reserve(seq.size());
    for(int x: seq){
        auto it = lower_bound(tail.begin(), tail.end(), x);
        if(it==tail.end()) tail.push_back(x);
        else *it = x;
    }
    return (int)tail.size();
}

// -------------------- Ulam distance via LIS on mapped positions --------------------
static int ulam_distance(const Perm& a, const Perm& b){
    // dU(a,b) = n - LCS(a,b)
    // For permutations, LCS(a,b) = LIS( pos_b[a[i]] )
    int n = (int)a.size();
    vector<int> pos(n+1);
    for(int i=0;i<n;i++) pos[b[i]] = i+1; // 1..n
    vector<int> seq(n);
    for(int i=0;i<n;i++) seq[i] = pos[a[i]];
    int L = lis_len_positions(seq);
    return n - L;
}

static bool is_valid_perm_1_to_n(const Perm& p){
    int n=(int)p.size();
    vector<int> seen(n+1,0);
    for(int x: p){
        if(x<1 || x>n) return false;
        if(seen[x]) return false;
        seen[x]=1;
    }
    return true;
}

// -------------------- IO helpers --------------------
static vector<Perm> read_perm_file(const string& path, int n){
    ifstream in(path);
    if(!in) throw runtime_error("Cannot open file: "+path);
    vector<Perm> perms;
    while(true){
        Perm p(n);
        for(int i=0;i<n;i++){
            if(!(in>>p[i])) return perms;
        }
        if(!is_valid_perm_1_to_n(p))
            throw runtime_error("Invalid permutation in file: "+path);
        perms.push_back(std::move(p));
    }
}

static void write_perm_file(const string& path, const vector<Perm>& perms){
    ofstream out(path);
    if(!out) throw runtime_error("Cannot write file: "+path);
    for(const auto& p: perms){
        for(size_t i=0;i<p.size();i++){
            out << p[i] << (i+1<p.size() ? ' ' : '\n');
        }
    }
}

static void write_one_perm(const string& path, const Perm& p){
    ofstream out(path);
    if(!out) throw runtime_error("Cannot write file: "+path);
    for(size_t i=0;i<p.size();i++){
        out << p[i] << (i+1<p.size() ? ' ' : '\n');
    }
}

static string perm_to_string(const Perm& p){
    string s;
    for(size_t i=0;i<p.size();i++){
        s += to_string(p[i]);
        if(i+1<p.size()) s += " ";
    }
    return s;
}

// -------------------- Build exact ball around identity using RSK+InverseRSK --------------------
static vector<Perm> build_ball_id(int n, int d){
    int r = d-1;
    int needFirst = n - r; // lambda_1 >= n-r

    auto shapes = gen_partitions(n, needFirst);

    vector<Perm> ball;
    // reserve rough amount if you know it; for safety:
    ball.reserve(5000);

    for(const auto& lam : shapes){
        // (optional) check hook-length count consistency
        cpp_int f = f_lambda(lam);
        auto syts = gen_SYT_from_shape(lam);
        if((cpp_int)syts.size() != f){
            // Not fatal, but should match
            // throw runtime_error("SYT count mismatch for a shape.");
        }

        for(const auto& P : syts){
            for(const auto& Q : syts){
                auto perm = inverse_rsk(P,Q);
                // perm is guaranteed in ball by construction
                ball.push_back(std::move(perm));
            }
        }
    }
    return ball;
}

// -------------------- Check candidate vs PA (distance>=d to all) --------------------
static bool far_from_all(const Perm& x, const vector<Perm>& PA, int d){
    for(const auto& s : PA){
        if(ulam_distance(x, s) < d) return false;
    }
    return true;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if(argc < 4){
        cerr << "Usage:\n";
        cerr << "  " << argv[0] << " n d pa_file [--max-trials T] [--seed S]\n\n";
        cerr << "Example:\n";
        cerr << "  " << argv[0] << " 11 3 11_3_pa.txt --max-trials 5000 --seed 1\n";
        return 1;
    }

    int n = stoi(argv[1]);
    int d = stoi(argv[2]);
    string pa_file = argv[3];

    long long max_trials = 2000;
    unsigned seed = (unsigned)chrono::high_resolution_clock::now().time_since_epoch().count();

    for(int i=4;i<argc;i++){
        string s = argv[i];
        if(s=="--max-trials" && i+1<argc) max_trials = stoll(argv[++i]);
        else if(s=="--seed" && i+1<argc) seed = (unsigned)stoul(argv[++i]);
    }

    // Read PA
    vector<Perm> PA = read_perm_file(pa_file, n);
    if(PA.empty()){
        cerr << "ERROR: PA file is empty.\n";
        return 1;
    }

    cout << "PA loaded: size = " << PA.size() << "\n";
    cout << "n="<<n<<", d="<<d<<", r=d-1="<<(d-1)<<"\n";

    // Build identity ball once
    cout << "Building B(id, r=d-1) using RSK+InverseRSK...\n";
    vector<Perm> ball_id = build_ball_id(n, d);
    cout << "|B(id, d-1)| = " << ball_id.size() << "\n";

    mt19937 rng(seed);
    uniform_int_distribution<size_t> dist_idx(0, PA.size()-1);

    // Try repeatedly
    for(long long trial=1; trial<=max_trials; trial++){
        // choose random center index
        size_t center_idx = dist_idx(rng);
        Perm center = PA[center_idx];

        // Build PA_minus = PA without center
        vector<Perm> PA_minus;
        PA_minus.reserve(PA.size()-1);
        for(size_t i=0;i<PA.size();i++){
            if(i!=center_idx) PA_minus.push_back(PA[i]);
        }

        // Collect candidates from ball(center)
        vector<Perm> candidates;
        candidates.reserve(16);

        for(const auto& tau : ball_id){
            Perm x = compose_left(center, tau);
            if(far_from_all(x, PA_minus, d)){
                candidates.push_back(std::move(x));
                if(candidates.size() >= 50) {
                    // keep bounded; we only need 2 with mutual dist>=d
                    // remove this cap if you want all candidates
                    break;
                }
            }
        }

        cout << "[trial " << trial << "] removed center, candidates found = " << candidates.size() << "\n";

        if(candidates.size() < 2) continue;

        // Need two candidates mutually distance >= d
        int a=-1,b=-1;
        for(int i=0;i<(int)candidates.size();i++){
            for(int j=i+1;j<(int)candidates.size();j++){
                if(ulam_distance(candidates[i], candidates[j]) >= d){
                    a=i; b=j; break;
                }
            }
            if(a!=-1) break;
        }

        if(a==-1){
            cout << "  candidates exist but no pair has mutual distance >= d\n";
            continue;
        }

        // SUCCESS: perform augmentation (remove 1 add 2)
        Perm add1 = candidates[a];
        Perm add2 = candidates[b];

        vector<Perm> PA_new = PA_minus;
        PA_new.push_back(add1);
        PA_new.push_back(add2);

        cout << "SUCCESS on trial " << trial << "!\n";
        cout << "Removed center: " << perm_to_string(center) << "\n";
        cout << "Add1:          " << perm_to_string(add1) << "\n";
        cout << "Add2:          " << perm_to_string(add2) << "\n";
        cout << "Old size="<<PA.size()<<", New size="<<PA_new.size()<<"  (improved by +1)\n";

        // Write outputs
        write_perm_file(pa_file, PA_new);
        write_one_perm("removed_center.txt", center);
        write_perm_file("added_two.txt", vector<Perm>{add1, add2});

        cout << "Updated PA written back to: " << pa_file << "\n";
        cout << "Saved removed center to: removed_center.txt\n";
        cout << "Saved added two perms to: added_two.txt\n";
        return 0;
    }

    cout << "No augmentation found within max_trials="<<max_trials<<". No files changed.\n";
    return 0;
}
