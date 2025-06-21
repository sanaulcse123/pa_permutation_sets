#include <bits/stdc++.h>
using namespace std;

vector<set<int>> graph;
vector<int> best_clique;

void bronKerbosch(set<int> R, set<int> P, set<int> X){
    if(P.empty() && X.empty()){
        if(R.size()>best_clique.size())
            best_clique=vector<int>(R.begin(),R.end());
        return;
    }
    int u=*P.begin();
    set<int> diff;
    for(int v:P)if(!graph[u].count(v))diff.insert(v);

    for(int v:diff){
        set<int>newR=R,newP,newX;
        newR.insert(v);
        for(int n:graph[v])if(P.count(n))newP.insert(n);
        for(int n:graph[v])if(X.count(n))newX.insert(n);
        bronKerbosch(newR,newP,newX);
        P.erase(v);
        X.insert(v);
    }
}

void read_dimacs(string fname,int &V){
    ifstream f(fname);
    string s;int u,v;
    while(f>>s)if(s=="p")break;
    f>>s>>V>>u;
    graph.resize(V+1);
    while(f>>s>>u>>v){
        graph[u].insert(v);graph[v].insert(u);
    }
    f.close();
}

void save_clique(string fname){
    ofstream f(fname);
    for(int v:best_clique)f<<v<<" ";
    f<<"\n";f.close();
}

void save_clique_perms(string permfile,string cliquefile,string fname){
    ifstream permf(permfile);
    vector<string> perms; string line;
    while(getline(permf,line)) perms.push_back(line);
    permf.close();

    ofstream outf(fname);
    for(int idx:best_clique) outf<<perms[idx-1]<<"\n";
    outf.close();
}

int main(){
    int n,d,V;
    cout<<"Enter n and d: ";
    cin>>n>>d;
    string dimacs="graph_n"+to_string(n)+"_d"+to_string(d)+".dimacs";
    string permfile="permutations_n"+to_string(n)+".txt";
    string cliquefile="clique_n"+to_string(n)+"_d"+to_string(d)+".clq";
    string permout="clique_n"+to_string(n)+"_d"+to_string(d)+".perm";

    read_dimacs(dimacs,V);
    set<int>P,X,R;
    for(int i=1;i<=V;i++)P.insert(i);
    bronKerbosch(R,P,X);
    
    save_clique(cliquefile);
    save_clique_perms(permfile,cliquefile,permout);

    cout<<" Clique size: "<<best_clique.size()<<endl;
}
