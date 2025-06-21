#include <bits/stdc++.h>
using namespace std;

int lcs(vector<int> &a, vector<int> &b) {
    int n = a.size();
    vector<vector<int>> dp(n+1, vector<int>(n+1,0));
    for(int i=1;i<=n;i++)
        for(int j=1;j<=n;j++)
            dp[i][j] = (a[i-1]==b[j-1]) ? dp[i-1][j-1]+1 : max(dp[i-1][j], dp[i][j-1]);
    return dp[n][n];
}

int levenshtein(vector<int>&a,vector<int>&b){
    return 2*(a.size()-lcs(a,b));
}

vector<vector<int>> generate_perms(int n){
    vector<int> p(n);
    iota(p.begin(),p.end(),1);
    vector<vector<int>> perms;
    do{ perms.push_back(p); }while(next_permutation(p.begin(),p.end()));
    return perms;
}

void save_perms(string fname,vector<vector<int>>&perms){
    ofstream f(fname);
    for(auto &p: perms){
        for(auto val:p)f<<val<<" ";
        f<<"\n";
    }
    f.close();
}

void save_dimacs(string fname,vector<pair<int,int>>& edges,int v){
    ofstream f(fname);
    f<<"p edge "<<v<<" "<<edges.size()<<"\n";
    for(auto e: edges)f<<"e "<<e.first<<" "<<e.second<<"\n";
    f.close();
}

int main(){
    int n,d;
    cout<<"Enter n and d: ";
    cin>>n>>d;
    auto perms=generate_perms(n);
    save_perms("permutations_n"+to_string(n)+".txt",perms);

    vector<pair<int,int>>edges;
    int V=perms.size();
    for(int i=0;i<V;i++){
        for(int j=i+1;j<V;j++){
            if(levenshtein(perms[i],perms[j])>=d)
                edges.push_back({i+1,j+1});
        }
    }
    save_dimacs("graph_n"+to_string(n)+"_d"+to_string(d)+".dimacs",edges,V);
    cout<<" Graph & permutation files created successfully!\n";
}
