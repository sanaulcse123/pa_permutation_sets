#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <fstream>

using namespace std;

// Efficient LCS computation
int lcs(const vector<int>& a, const vector<int>& b) {
    int n = a.size();
    vector<int> dp(n + 1, 0), prev(n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        swap(dp, prev);
        for (int j = 1; j <= n; ++j) {
            dp[j] = (a[i - 1] == b[j - 1]) ? prev[j - 1] + 1 : max(prev[j], dp[j - 1]);
        }
    }
    return dp[n];
}

// Levenshtein = 2 * (n - LCS)
int levenshtein(const vector<int>& a, const vector<int>& b) {
    return 2 * (a.size() - lcs(a, b));
}

// Generate random permutation
vector<int> generatePermutation(int n, mt19937& rng) {
    vector<int> perm(n);
    for (int i = 0; i < n; ++i) perm[i] = i + 1;
    shuffle(perm.begin(), perm.end(), rng);
    return perm;
}

int main() {
    int n, d;
    cout << "Enter n : ";
    cin >> n;
    cout << "Enter d : ";
    cin >> d;

    const int TRIALS = 500000;
    const int OUTER = 100;
    vector<vector<int>> bestSet;
    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

    for (int attempt = 0; attempt < OUTER; ++attempt) {
        vector<vector<int>> currentSet;

        for (int i = 0; i < TRIALS; ++i) {
            vector<int> candidate = generatePermutation(n, rng);
            bool valid = true;

            for (const auto& p : currentSet) {
                if (levenshtein(candidate, p) < d) {
                    valid = false;
                    break;
                }
            }

            if (valid)
                currentSet.push_back(candidate);
        }

        if (currentSet.size() > bestSet.size()) {
            bestSet = currentSet;
            cout << "[!] New best found in attempt " << attempt + 1
                 << ": Size = " << bestSet.size() << endl;
        }
    }

    cout << "\n Estimated L(" << n << ", " << d << ") = " << bestSet.size() << endl;

    // Save result
    ofstream fout("best_L_" + to_string(n) + "_" + to_string(d) + ".txt");
    for (const auto& p : bestSet) {
        for (int x : p) fout << x << " ";
        fout << "\n";
    }
    fout.close();

    return 0;
}
