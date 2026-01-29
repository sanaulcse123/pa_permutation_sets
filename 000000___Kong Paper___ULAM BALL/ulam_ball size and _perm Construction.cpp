// 5_3perm.cpp  (Paper-style construction)
// Input: n d
// Output: For all partitions (shapes) λ of n with λ1 >= n-(d-1),
// generate ALL permutations of shape λ using:
//    all SYT P of shape λ  ×  all SYT Q of shape λ  → inverse RSK → permutation
// Hence exactly (f^λ)^2 permutations per λ.
//
// Works great for (n,d)=(5,3): shapes (5),(4,1),(3,2),(3,1,1) with counts 1,16,25,36.

#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using boost::multiprecision::cpp_int;

// ---------- helpers ----------
static cpp_int factorial_int(int n){
    cpp_int f=1;
    for(int i=2;i<=n;i++) f*=i;
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
    // sort by first row desc, then lex desc (nice printing)
    sort(out.begin(), out.end(), [](const auto& a, const auto& b){
        if(a[0]!=b[0]) return a[0]>b[0];
        return a>b;
    });
    return out;
}

// ---------- SYT generation by removing corners (largest number must be in a corner) ----------
// Represent tableau as vector< vector<int> > rows, each row length = shape row length.

using Tableau = vector<vector<int>>;

static vector<int> corners_of_shape(const vector<int>& lam){
    // corner row indices i where cell (i, lam[i]-1) is a corner:
    // i is corner if i is last row OR lam[i] > lam[i+1]
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
    // returns all SYT of shape lam filled with numbers 1..k
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
            // add back the removed corner cell and place k there
            if(rowRemoved){
                // insert a new row at position cr with one cell = k
                T.insert(T.begin()+cr, vector<int>{k});
            }else{
                // append to end of existing row cr
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
// Given SYT P,Q (same shape) produce permutation w[1..n].
static vector<int> inverse_rsk(Tableau P, Tableau Q){
    int n=0;
    for(auto &row: Q) n += (int)row.size();

    vector<int> w(n+1, 0); // 1-indexed

    for(int k=n; k>=1; --k){
        // find position of k in Q
        int r=-1, c=-1;
        for(int i=0;i<(int)Q.size();i++){
            for(int j=0;j<(int)Q[i].size();j++){
                if(Q[i][j]==k){ r=i; c=j; break; }
            }
            if(r!=-1) break;
        }
        if(r==-1){
            throw runtime_error("inverse_rsk: cannot find k in Q");
        }

        // remove cell (r,c) from Q (it should be a corner during correct inversion)
        Q[r].erase(Q[r].begin()+c);
        if(Q[r].empty()) Q.erase(Q.begin()+r);

        // remove same cell from P and take its value x
        int x = P[r][c];
        P[r].erase(P[r].begin()+c);
        if(P[r].empty()) P.erase(P.begin()+r);

        // reverse bump upwards: for rows r-1 down to 0
        for(int i=r-1; i>=0; --i){
            auto &row = P[i];
            // find rightmost element < x
            auto it = lower_bound(row.begin(), row.end(), x); // first >= x
            if(it==row.begin()){
                throw runtime_error("inverse_rsk: reverse bump failed (no element < x)");
            }
            --it;
            int y = *it;
            *it = x;
            x = y;
        }
        w[k] = x;
    }

    // convert to 0-indexed vector<int> of size n
    vector<int> perm(n);
    for(int i=1;i<=n;i++) perm[i-1]=w[i];
    return perm;
}

// ---------- Ulam distance to id via LIS (for verification only) ----------
static int lis_len(const vector<int>& p){
    vector<int> tail;
    for(int x: p){
        auto it = lower_bound(tail.begin(), tail.end(), x);
        if(it==tail.end()) tail.push_back(x);
        else *it = x;
    }
    return (int)tail.size();
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n,d;
    if(argc>=3){
        n = stoi(argv[1]);
        d = stoi(argv[2]);
    }else{
        cout << "Enter n d (example: 5 3): ";
        cin >> n >> d;
    }

    int r = d-1;
    int needFirst = n - r; // λ1 >= needFirst

    cout << "\n=== PAPER-style construction: (P,Q) SYT pairs -> inverse RSK -> permutations ===\n";
    cout << "n="<<n<<", d="<<d<<", radius r=d-1="<<r<<", need LIS = needFirst = "<<needFirst<<"\n\n";

    auto shapes = gen_partitions(n, needFirst);

    long long grandTotal = 0;

    for(const auto& lam : shapes){
        cpp_int f = f_lambda(lam);
        cpp_int f2 = f*f;

        auto syts = gen_SYT_from_shape(lam);

        cout << "------------------------------------------\n";
        cout << "Shape lambda = " << shape_to_string(lam) << "\n";
        cout << "f^lambda = " << f << "   =>   (f^lambda)^2 = " << f2 << "\n";
        cout << "SYT count generated = " << syts.size() << "\n";

        // Generate permutations from all pairs (P,Q)
        vector<vector<int>> perms;
        perms.reserve((size_t)syts.size() * (size_t)syts.size());

        for(const auto& P : syts){
            for(const auto& Q : syts){
                auto perm = inverse_rsk(P,Q);
                // optional verification: dist <= r
                int L = lis_len(perm);
                int dist = n - L;
                if(dist > r){
                    cerr << "ERROR: generated perm not in ball? dist="<<dist<<"\n";
                }
                perms.push_back(std::move(perm));
            }
        }

        cout << "Permutations count = " << perms.size() << "\n";
        cout << "Permutations:\n";
        for(const auto& perm : perms){
            for(int i=0;i<n;i++){
                cout << perm[i] << (i+1<n ? ' ' : '\n');
            }
        }
        cout << "\n";

        grandTotal += (long long)perms.size();
    }

    cout << "Total permutations in ball B("<<n<<","<<r<<") = " << grandTotal << "\n";
    cout << "NOTE: This is exact paper-style construction. It explodes fast for large n.\n";
    return 0;
}
