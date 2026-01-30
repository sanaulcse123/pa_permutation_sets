// n_d_ball.cpp
// Exact construction of Ulam ball B(id, r=d-1) using:
//   shapes λ with λ1 >= n-r, all SYT pairs (P,Q), InverseRSK(P,Q).
// Optional center shift using left-translation:
//   B(pi,r) = { pi ∘ tau : tau ∈ B(id,r) }.
//
// Build (WSL):
//   g++ -O2 -std=c++17 n_d_ball.cpp -o n_d_ball
//
// Run:
//   ./n_d_ball 5 3
//   ./n_d_ball 5 3 --center 1 3 5 2 4

#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using boost::multiprecision::cpp_int;

using Tableau = vector<vector<int>>;

// -------------------- small helpers --------------------
static cpp_int factorial_int(int n){
    cpp_int f = 1;
    for(int i=2;i<=n;i++) f *= i;
    return f;
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

static void print_perm(const vector<int>& p){
    for(size_t i=0;i<p.size();i++){
        cout << p[i] << (i+1<p.size() ? ' ' : '\n');
    }
}

// σ = π ∘ τ in one-line notation: σ_i = π_{τ_i}
static vector<int> compose_left(const vector<int>& pi, const vector<int>& tau){
    int n = (int)pi.size();
    vector<int> sigma(n);
    for(int i=0;i<n;i++) sigma[i] = pi[tau[i]-1];
    return sigma;
}

// -------------------- Hook-length f^λ (for reporting only) --------------------
static cpp_int f_lambda(const vector<int>& lam){
    int n=0; for(int x: lam) n+=x;
    cpp_int num = factorial_int(n);
    cpp_int den = 1;

    for(int i=0;i<(int)lam.size();i++){
        for(int j=0;j<lam[i];j++){
            int right = lam[i] - (j+1);
            int below = 0;
            for(int k=i+1;k<(int)lam.size();k++){
                if(lam[k] > j) below++;
            }
            den *= (right + below + 1);
        }
    }
    return num/den;
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

// -------------------- generate partitions λ ⊢ n with λ1 >= minFirst --------------------
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

// -------------------- AllSYT(λ) by removing corners --------------------
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

// -------------------- InverseRSK(P,Q) (reverse row insertion) --------------------
static vector<int> InverseRSK(Tableau P, Tableau Q){
    int n=0; for(auto &row: Q) n += (int)row.size();
    vector<int> w(n+1, 0); // 1-indexed

    for(int k=n; k>=1; --k){
        int r=-1, c=-1;
        for(int i=0;i<(int)Q.size();i++){
            for(int j=0;j<(int)Q[i].size();j++){
                if(Q[i][j]==k){ r=i; c=j; break; }
            }
            if(r!=-1) break;
        }
        if(r==-1) throw runtime_error("InverseRSK: cannot find k in Q");

        // remove that cell from Q
        Q[r].erase(Q[r].begin()+c);
        if(Q[r].empty()) Q.erase(Q.begin()+r);

        // remove same cell from P and reverse-bump upward
        int x = P[r][c];
        P[r].erase(P[r].begin()+c);
        if(P[r].empty()) P.erase(P.begin()+r);

        for(int i=r-1; i>=0; --i){
            auto &row = P[i];
            auto it = lower_bound(row.begin(), row.end(), x); // first >= x
            if(it==row.begin()) throw runtime_error("InverseRSK: reverse bump failed");
            --it;               // rightmost < x
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

// -------------------- BallConstruct (Section 6) --------------------
static vector<vector<int>> BallConstruct(int n, int d, bool print_by_shape=true){
    int r = d-1;
    int minFirst = n - r; // λ1 >= n-r
    auto shapes = gen_partitions(n, minFirst);

    vector<vector<int>> ball_id;
    cpp_int total = 0;

    cout << "n="<<n<<", d="<<d<<", r=d-1="<<r<<", condition: LIS >= "<<minFirst<<"\n\n";

    for(const auto& lam : shapes){
        auto T = AllSYT(lam);
        cpp_int f = f_lambda(lam);

        cout << "------------------------------------------\n";
        cout << "Shape lambda = " << shape_to_string(lam) << "\n";
        cout << "f^lambda = " << f << "   =>   (f^lambda)^2 = " << (f*f) << "\n";
        cout << "SYT count = " << T.size() << "\n";

        vector<vector<int>> perms_shape;
        perms_shape.reserve((size_t)T.size() * (size_t)T.size());

        for(const auto& P : T){
            for(const auto& Q : T){
                perms_shape.push_back(InverseRSK(P,Q));
            }
        }

        cout << "Permutations count = " << perms_shape.size() << "\n";
        if(print_by_shape){
            cout << "Permutations:\n";
            for(const auto& p : perms_shape) print_perm(p);
            cout << "\n";
        }

        total += (cpp_int)perms_shape.size();
        // append to global ball list
        for(auto &p : perms_shape) ball_id.push_back(std::move(p));
    }

    cout << "Total |B(id,r)| = " << total << "\n\n";
    return ball_id;
}

// -------------------- TranslateBall (Section 6) --------------------
static void TranslateBallPrint(const vector<int>& pi, const vector<vector<int>>& ball_id){
    cout << "==========================================\n";
    cout << "Ball centered at pi:  B(pi,r) = { pi ∘ tau : tau in B(id,r) }\n";
    cout << "pi = "; print_perm(pi);
    cout << "Ball size = " << ball_id.size() << "\n";
    cout << "Permutations in B(pi,r):\n";
    for(const auto& tau : ball_id){
        auto sigma = compose_left(pi, tau);
        print_perm(sigma);
    }
    cout << "==========================================\n";
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n=-1, d=-1;
    bool have_center=false;
    vector<int> center_pi;

    if(argc >= 3){
        n = stoi(argv[1]);
        d = stoi(argv[2]);
        for(int i=3;i<argc;i++){
            if(string(argv[i])=="--center"){
                have_center=true;
                center_pi.clear();
                for(int j=i+1;j<argc;j++) center_pi.push_back(stoi(argv[j]));
                break;
            }
        }
    }else{
        cout << "Enter n d (example: 5 3): ";
        cin >> n >> d;
        cout << "Compute ball centered at another permutation? (y/n): ";
        char ch; cin >> ch;
        if(ch=='y' || ch=='Y'){
            have_center=true;
            center_pi.assign(n,0);
            cout << "Enter center permutation pi ("<<n<<" numbers): ";
            for(int i=0;i<n;i++) cin >> center_pi[i];
        }
    }

    if(n<=0 || d<=0){
        cerr << "ERROR: invalid n or d\n";
        return 1;
    }
    if(have_center){
        if((int)center_pi.size()!=n){
            cerr << "ERROR: --center needs exactly n numbers\n";
            return 1;
        }
        if(!is_valid_perm_1_to_n(center_pi)){
            cerr << "ERROR: center is not a permutation of 1..n\n";
            return 1;
        }
    }

    cout << "=== PAPER-style construction: (P,Q) SYT pairs -> InverseRSK -> permutations ===\n";
    auto ball_id = BallConstruct(n,d,true);

    if(have_center){
        TranslateBallPrint(center_pi, ball_id);
    }

    return 0;
}
