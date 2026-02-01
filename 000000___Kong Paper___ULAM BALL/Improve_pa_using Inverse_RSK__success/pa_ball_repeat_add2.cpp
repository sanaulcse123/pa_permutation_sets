// pa_ball_repeat_add2_safe.cpp
// Remove 1 permutation pi from an (n,d)-PA under Ulam metric,
// search inside the Ulam ball B(pi, r=d-1) for TWO permutations add1,add2 such that:
//   (1) dU(add1, sigma) >= d for all sigma in PA \ {pi}
//   (2) dU(add2, sigma) >= d for all sigma in PA \ {pi}
//   (3) dU(add1, add2)  >= d
// If found: output improved PA of size |PA|+1 (remove 1, add 2).
//
// IMPORTANT: We do NOT require distance >= d to the removed center pi (impossible since add's are in B(pi,d-1)).
//
// Ball construction: exact B(id, d-1) via (P,Q) SYT pairs + inverse RSK, then translate: B(pi,r)= { pi ∘ tau }.
//
// SAFE behavior: PA file is written ONLY on success. Otherwise unchanged.
// Also makes a backup <pa_file>.bak on success.

#include <bits/stdc++.h>
using namespace std;

using Perm = vector<int>;
using Tableau = vector<vector<int>>;

// ------------------ small helpers ------------------
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

static string perm_to_string(const Perm& p){
    string s;
    for(size_t i=0;i<p.size();i++){
        s += to_string(p[i]);
        if(i+1<p.size()) s += " ";
    }
    return s;
}

static void write_perm_file(const string& path, const vector<Perm>& perms){
    ofstream out(path);
    if(!out) throw runtime_error("Cannot write file: " + path);
    for(const auto& p: perms){
        for(size_t i=0;i<p.size();i++){
            out << p[i] << (i+1<p.size() ? ' ' : '\n');
        }
    }
}

static void write_one_perm(const string& path, const Perm& p){
    ofstream out(path);
    if(!out) throw runtime_error("Cannot write file: " + path);
    for(size_t i=0;i<p.size();i++){
        out << p[i] << (i+1<p.size() ? ' ' : '\n');
    }
}

static vector<Perm> read_perm_file(const string& path, int n){
    ifstream in(path);
    if(!in) throw runtime_error("Cannot open file: " + path);

    vector<Perm> perms;
    while(true){
        Perm p(n);
        for(int i=0;i<n;i++){
            if(!(in >> p[i])) return perms; // normal EOF
        }
        if(!is_valid_perm_1_to_n(p))
            throw runtime_error("Invalid permutation found in: " + path);
        perms.push_back(std::move(p));
    }
}

// ------------------ partitions of n with lambda1 >= minFirst ------------------
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

// ------------------ SYT generation by removing corners ------------------
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

// ------------------ Inverse RSK (reverse row insertion) ------------------
static Perm inverse_rsk(Tableau P, Tableau Q){
    int n=0;
    for(auto &row: Q) n += (int)row.size();

    vector<int> w(n+1, 0); // 1-indexed output

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

        // reverse bump up
        for(int i=r-1; i>=0; --i){
            auto &row = P[i];
            auto it = lower_bound(row.begin(), row.end(), x); // first >= x
            if(it==row.begin()) throw runtime_error("inverse_rsk: reverse bump failed");
            --it;                 // rightmost < x
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

// ------------------ composition (left translate): pi ∘ tau ------------------
static Perm compose_left(const Perm& pi, const Perm& tau){
    int n = (int)pi.size();
    Perm sigma(n);
    for(int i=0;i<n;i++){
        sigma[i] = pi[tau[i]-1];
    }
    return sigma;
}

// ------------------ Ulam distance fast: dU(a,b)=n-LIS(pos_b[a[i]]) ------------------
static int lis_len_n2(const vector<int>& seq){
    int n = (int)seq.size();
    vector<int> dp(n,1);
    int best=0;
    for(int i=0;i<n;i++){
        for(int j=0;j<i;j++){
            if(seq[j] < seq[i]) dp[i] = max(dp[i], dp[j]+1);
        }
        best = max(best, dp[i]);
    }
    return best;
}

static vector<uint16_t> build_pos(const Perm& b){
    int n = (int)b.size();
    vector<uint16_t> pos(n+1);
    for(int i=0;i<n;i++) pos[b[i]] = (uint16_t)(i+1); // 1..n
    return pos;
}

static int ulam_distance_with_pos(const Perm& a, const vector<uint16_t>& pos_b){
    int n = (int)a.size();
    vector<int> seq(n);
    for(int i=0;i<n;i++) seq[i] = (int)pos_b[a[i]];
    int L = lis_len_n2(seq);
    return n - L;
}

// ------------------ Build exact B(id, d-1) once ------------------
static vector<Perm> build_ball_id(int n, int d){
    int r = d-1;
    int needFirst = n - r; // lambda1 >= n-r
    auto shapes = gen_partitions(n, needFirst);

    vector<Perm> ball;
    // for (11,3): size is 4062, so safe reserve:
    ball.reserve(8000);

    for(const auto& lam : shapes){
        auto syts = gen_SYT_from_shape(lam);
        for(const auto& P : syts){
            for(const auto& Q : syts){
                ball.push_back(inverse_rsk(P,Q));
            }
        }
    }

    // optional dedup (should already be unique):
    // sort and unique by string key for safety
    vector<string> keys;
    keys.reserve(ball.size());
    for(auto &p: ball) keys.push_back(perm_to_string(p));
    vector<int> idx(ball.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int i, int j){ return keys[i] < keys[j]; });

    vector<Perm> uniq;
    uniq.reserve(ball.size());
    for(size_t t=0;t<idx.size();t++){
        if(t==0 || keys[idx[t]] != keys[idx[t-1]]) uniq.push_back(ball[idx[t]]);
    }
    return uniq;
}

// check x far from all PA entries except center_idx
static bool far_from_all_except(const Perm& x,
                                const vector<Perm>& PA,
                                const vector<vector<uint16_t>>& PApos,
                                size_t center_idx,
                                int d){
    for(size_t i=0;i<PA.size();i++){
        if(i==center_idx) continue;
        if(ulam_distance_with_pos(x, PApos[i]) < d) return false;
    }
    return true;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if(argc < 4){
        cerr << "Usage:\n";
        cerr << "  " << argv[0] << " n d pa_file [--out OUTFILE] [--max-trials T] [--seed S] [--candidate-limit K] [--print-every P]\n\n";
        cerr << "Example:\n";
        cerr << "  " << argv[0] << " 11 3 11_3_pa.txt --max-trials 20000 --seed 1 --candidate-limit 500 --out 11_3_pa_improved.txt\n";
        return 1;
    }

    int n = stoi(argv[1]);
    int d = stoi(argv[2]);
    string pa_file = argv[3];

    long long max_trials = 20000;
    unsigned seed = (unsigned)chrono::high_resolution_clock::now().time_since_epoch().count();
    size_t candidate_limit = 500;   // collect up to K candidates per trial (increase if you want)
    long long print_every = 1;       // print every trial
    string out_file = pa_file;       // default overwrite on success

    for(int i=4;i<argc;i++){
        string s = argv[i];
        if(s=="--max-trials" && i+1<argc) max_trials = stoll(argv[++i]);
        else if(s=="--seed" && i+1<argc) seed = (unsigned)stoul(argv[++i]);
        else if(s=="--candidate-limit" && i+1<argc) candidate_limit = (size_t)stoull(argv[++i]);
        else if(s=="--print-every" && i+1<argc) print_every = stoll(argv[++i]);
        else if(s=="--out" && i+1<argc) out_file = argv[++i];
    }

    // 1) read PA
    vector<Perm> PA = read_perm_file(pa_file, n);
    if(PA.empty()){
        cerr << "ERROR: PA file is empty or not readable as integers.\n";
        cerr << "Tip: run `wc -l " << pa_file << "` and `head -n 2 " << pa_file << "`\n";
        return 1;
    }

    // precompute pos arrays for PA
    vector<vector<uint16_t>> PApos(PA.size());
    for(size_t i=0;i<PA.size();i++) PApos[i] = build_pos(PA[i]);

    cout << "PA loaded: size = " << PA.size() << "\n";
    cout << "n="<<n<<", d="<<d<<", r=d-1="<<(d-1)<<"\n";
    cout << "Seed="<<seed<<", max_trials="<<max_trials<<", candidate_limit="<<candidate_limit<<"\n";

    // 2) build B(id, d-1)
    cout << "Building exact B(id, r=d-1) using SYT pairs + inverse RSK...\n";
    vector<Perm> ball_id = build_ball_id(n, d);
    cout << "|B(id, d-1)| = " << ball_id.size() << "\n";

    mt19937 rng(seed);
    uniform_int_distribution<size_t> dist_idx(0, PA.size()-1);

    // 3) trials (NO file write until success)
    for(long long trial=1; trial<=max_trials; trial++){
        size_t center_idx = dist_idx(rng);
        const Perm& center = PA[center_idx];

        // collect candidates in B(center)
        vector<Perm> candidates;
        candidates.reserve(min(candidate_limit, (size_t)ball_id.size()));

        for(const auto& tau : ball_id){
            Perm x = compose_left(center, tau); // x in B(center)
            if(far_from_all_except(x, PA, PApos, center_idx, d)){
                candidates.push_back(std::move(x));
                if(candidates.size() >= candidate_limit) break;
            }
        }

        if(print_every>0 && (trial % print_every)==0){
            cout << "[trial " << trial << "] candidates found = " << candidates.size() << "\n";
        }

        if(candidates.size() < 2) continue;

        // need two candidates mutually distance >= d
        vector<vector<uint16_t>> candPos(candidates.size());
        for(size_t i=0;i<candidates.size();i++) candPos[i] = build_pos(candidates[i]);

        int a=-1,b=-1;
        for(int i=0;i<(int)candidates.size();i++){
            for(int j=i+1;j<(int)candidates.size();j++){
                int dist_ab = ulam_distance_with_pos(candidates[i], candPos[j]);
                if(dist_ab >= d){
                    a=i; b=j; break;
                }
            }
            if(a!=-1) break;
        }
        if(a==-1) continue;

        // FINAL correctness re-check before writing:
        // add1 far from all except center, add2 far from all except center, and mutual ok.
        Perm add1 = candidates[a];
        Perm add2 = candidates[b];

        if(!far_from_all_except(add1, PA, PApos, center_idx, d)) continue;
        if(!far_from_all_except(add2, PA, PApos, center_idx, d)) continue;

        auto add2pos = build_pos(add2);
        if(ulam_distance_with_pos(add1, add2pos) < d) continue;

        // SUCCESS -> build PA_new
        vector<Perm> PA_new;
        PA_new.reserve(PA.size() + 1);
        for(size_t i=0;i<PA.size();i++){
            if(i==center_idx) continue;
            PA_new.push_back(PA[i]);
        }
        PA_new.push_back(add1);
        PA_new.push_back(add2);

        cout << "SUCCESS on trial " << trial << "!\n";
        cout << "Removed center: " << perm_to_string(center) << "\n";
        cout << "Add1:          " << perm_to_string(add1) << "\n";
        cout << "Add2:          " << perm_to_string(add2) << "\n";
        cout << "Old size=" << PA.size() << ", New size=" << PA_new.size() << " (improved by +1)\n";

        // backup original pa_file
        {
            string bak = pa_file + ".bak";
            // overwrite backup each success (simple & safe)
            ifstream src(pa_file, ios::binary);
            ofstream dst(bak, ios::binary);
            dst << src.rdbuf();
            cout << "Backup saved: " << bak << "\n";
        }

        // write outputs
        write_perm_file(out_file, PA_new);
        write_one_perm("removed_center.txt", center);
        write_perm_file("added_two.txt", vector<Perm>{add1, add2});

        cout << "Updated PA written to: " << out_file << "\n";
        cout << "Saved removed center to: removed_center.txt\n";
        cout << "Saved added two perms to: added_two.txt\n";
        return 0;
    }

    cout << "No augmentation found within max_trials="<<max_trials<<". (No files changed.)\n";
    return 0;
}
