// n_d_ball.cpp
// Compute exact Ulam ball B(id, r=d-1) in S_n using (P,Q) SYT pairs -> inverse RSK,
// and (optionally) compute the ball centered at an arbitrary permutation pi by left-translation:
//      B(pi, r) = { pi ∘ tau : tau ∈ B(id, r) }.
//
// Usage (WSL):
//   g++ -O2 -std=c++17 n_d_ball.cpp -o n_d_ball
//   ./n_d_ball 5 3
//   ./n_d_ball 5 3 --center 1 3 5 2 4
//
// If you don't pass args, it will prompt interactively.

#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using boost::multiprecision::cpp_int;

// ---------- helpers ----------
static cpp_int factorial_int(int n){
    cpp_int f = 1;
    for(int i=2;i<=n;i++) f *= i;
    return f;
}

static string shape_to_string(const vector<int>& lam){
    string s="(";
    for(size_t i=0;i<lam.size();i++){
        s += to_string(lam[i]);
        if(i+1<lam.size()) s += ",";
    }
    s += ")";
    return s;
}

// Hook-length formula: f^lambda = n! / prod hooks
static cpp_int f_lambda(const vector<int>& lam){
    int n=0;
    for(int x: lam) n+=x;

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

// ---------- generate all partitions of n with first part >= minFirst ----------
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

// ---------- SYT generation by removing corners ----------
using Tableau = vector<vector<int>>;

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
    if(lam2[cornerRow]==0){
        lam2.erase(lam2.begin()+cornerRow);
    }
    return lam2;
}

static vector<Tableau> gen_SYT_from_shape_rec(const vector<int>& lam, int k){
    if(k==0){
        return vector<Tableau>{Tableau{}};
    }
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

// ---------- Inverse RSK (reverse row insertion) ----------
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

        for(int i=r-1; i>=0; --i){
            auto &row = P[i];
            auto it = lower_bound(row.begin(), row.end(), x); // first >= x
            if(it==row.begin()){
                throw runtime_error("inverse_rsk: reverse bump failed (no element < x)");
            }
            --it;                 // rightmost < x
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

// ---------- LIS length (verification) ----------
static int lis_len(const vector<int>& p){
    vector<int> tail;
    for(int x: p){
        auto it = lower_bound(tail.begin(), tail.end(), x);
        if(it==tail.end()) tail.push_back(x);
        else *it = x;
    }
    return (int)tail.size();
}

// ---------- Compose (left translate): sigma = pi ∘ tau ----------
static vector<int> compose_left(const vector<int>& pi, const vector<int>& tau){
    int n = (int)pi.size();
    vector<int> sigma(n);
    for(int i=0;i<n;i++){
        sigma[i] = pi[tau[i]-1]; // one-line notation
    }
    return sigma;
}

static void print_perm(const vector<int>& p){
    for(size_t i=0;i<p.size();i++){
        cout << p[i] << (i+1<p.size() ? ' ' : '\n');
    }
}

static bool is_valid_perm_1_to_n(const vector<int>& p){
    int n = (int)p.size();
    vector<int> seen(n+1,0);
    for(int x: p){
        if(x<1 || x>n) return false;
        if(seen[x]) return false;
        seen[x]=1;
    }
    return true;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n=-1, d=-1;
    vector<int> center_pi;
    bool have_center = false;

    // ------------------ parse args ------------------
    // Supported:
    //   ./n_d_ball n d
    //   ./n_d_ball n d --center p1 p2 ... pn
    if(argc >= 3){
        n = stoi(argv[1]);
        d = stoi(argv[2]);
        for(int i=3;i<argc;i++){
            string s = argv[i];
            if(s=="--center"){
                have_center = true;
                center_pi.clear();
                for(int j=i+1; j<argc; j++){
                    center_pi.push_back(stoi(argv[j]));
                }
                break;
            }
        }
        if(have_center && (int)center_pi.size() != n){
            cerr << "ERROR: --center provided but number of elements is " << center_pi.size()
                 << " (expected n=" << n << ")\n";
            return 1;
        }
    }else{
        cout << "Enter n d (example: 5 3): ";
        cin >> n >> d;
        cout << "Compute ball centered at another permutation? (y/n): ";
        char ch; cin >> ch;
        if(ch=='y' || ch=='Y'){
            have_center = true;
            center_pi.assign(n,0);
            cout << "Enter center permutation pi ("<<n<<" numbers): ";
            for(int i=0;i<n;i++) cin >> center_pi[i];
        }
    }

    if(n<=0 || d<=0){
        cerr << "ERROR: invalid n or d\n";
        return 1;
    }
    if(have_center && !is_valid_perm_1_to_n(center_pi)){
        cerr << "ERROR: center pi is not a valid permutation of 1..n\n";
        return 1;
    }

    int r = d-1;
    int needFirst = n - r; // λ1 >= needFirst

    cout << "\n=== PAPER-style construction: (P,Q) SYT pairs -> inverse RSK -> permutations ===\n";
    cout << "n="<<n<<", d="<<d<<", radius r=d-1="<<r<<", need LIS = "<<needFirst<<"\n\n";

    auto shapes = gen_partitions(n, needFirst);

    // We'll store the entire ball around identity here:
    vector<vector<int>> ball_id;
    ball_id.reserve(1000);

    cpp_int grandTotal = 0;

    for(const auto& lam : shapes){
        cpp_int f = f_lambda(lam);
        cpp_int f2 = f*f;

        auto syts = gen_SYT_from_shape(lam);

        cout << "------------------------------------------\n";
        cout << "Shape lambda = " << shape_to_string(lam) << "\n";
        cout << "f^lambda = " << f << "   =>   (f^lambda)^2 = " << f2 << "\n";
        cout << "SYT count generated = " << syts.size() << "\n";

        cpp_int count_here = 0;

        // Generate permutations from all pairs (P,Q)
        for(const auto& P : syts){
            for(const auto& Q : syts){
                auto perm = inverse_rsk(P,Q);

                // verification (same as your existing code)
                int L = lis_len(perm);
                int dist = n - L;
                if(dist > r){
                    cerr << "ERROR: generated perm not in ball? dist="<<dist<<"\n";
                    return 1;
                }

                ball_id.push_back(std::move(perm));
                count_here += 1;
            }
        }

        cout << "Permutations count = " << count_here << "\n";
        cout << "Permutations:\n";
        // Print only this shape's block (last count_here perms)
        // (so output format stays like your existing program)
        size_t start = ball_id.size() - (size_t)count_here.convert_to<long long>();
        for(size_t idx=start; idx<ball_id.size(); idx++){
            print_perm(ball_id[idx]);
        }
        cout << "\n";

        grandTotal += count_here;
    }

    cout << "Total permutations in ball B("<<n<<","<<r<<") = " << grandTotal << "\n";
    cout << "NOTE: Exact paper-style construction (no brute force over n!).\n\n";

    // ---------- NEW FEATURE: ball centered at pi ----------
    if(have_center){
        cout << "==========================================\n";
        cout << "Ball centered at pi (same radius r="<<r<<")\n";
        cout << "pi = "; print_perm(center_pi);
        cout << "Construction: B(pi,r) = { pi ∘ tau : tau ∈ B(id,r) }\n";
        cout << "Ball size = " << ball_id.size() << "\n";
        cout << "Permutations in B(pi,r):\n";

        for(const auto& tau : ball_id){
            auto sigma = compose_left(center_pi, tau);
            print_perm(sigma);
        }
        cout << "==========================================\n";
    }

    return 0;
}
