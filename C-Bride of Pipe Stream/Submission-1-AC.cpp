#include <bits/stdc++.h>
using namespace std;

struct Duct {
    vector<pair<int,double>> toStation; // (station index, fraction)
    vector<double> toRes; // size r
};

int s, r, dnum;
vector<vector<Duct>> outs;

double evalH1(const vector<double>& lam){
    vector<double> H(s+1, 0.0);
    for(int i = s; i >= 1; --i){
        double best = -1e100;
        for(const Duct& duct : outs[i]){
            double val = 0.0;
            for(int j=0;j<r;++j) val += duct.toRes[j] * lam[j];
            for(auto &pr : duct.toStation){
                val += pr.second * H[pr.first];
            }
            if(val > best) best = val;
        }
        H[i] = best;
    }
    return H[1];
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    if(!(cin>>s>>r>>dnum)) return 0;
    outs.assign(s+1, {});
    for(int di=0; di<dnum; ++di){
        int i, n; cin>>i>>n;
        Duct duct; duct.toStation.clear(); duct.toRes.assign(r, 0.0);
        for(int k=0;k<n;++k){
            int o, p; cin>>o>>p; double frac = p/100.0;
            if(o <= s) duct.toStation.push_back({o, frac});
            else {
                int j = o - s - 1; 
                if(j>=0 && j<r) duct.toRes[j] += frac;
            }
        }
        outs[i].push_back(move(duct));
    }
    cout.setf(ios::fixed);
    cout<<setprecision(10);
    if(r==1){
        vector<double> lam = {1.0};
        double ans = evalH1(lam) * 100.0;
        cout<<ans<<"\n";
        return 0;
    } else if(r==2){
        double l=0.0, rgt=1.0;
        const double gr = (sqrt(5.0)-1.0)/2.0;
        double x1 = rgt - gr * (rgt - l);
        double x2 = l + gr * (rgt - l);
        vector<double> lam1={x1,1.0-x1}, lam2={x2,1.0-x2};
        double f1 = evalH1(lam1), f2 = evalH1(lam2);
        for(int it=0; it<200; ++it){
            if(f1 > f2){
                l = x1;
                x1 = x2;
                f1 = f2;
                x2 = l + gr*(rgt - l);
                lam2[0] = x2; lam2[1] = 1.0 - x2;
                f2 = evalH1(lam2);
            } else {
                rgt = x2;
                x2 = x1;
                f2 = f1;
                x1 = rgt - gr*(rgt - l);
                lam1[0] = x1; lam1[1] = 1.0 - x1;
                f1 = evalH1(lam1);
            }
        }
        double ans = min(f1, f2)*100.0;
        cout<<ans<<"\n";
        return 0;
    } else { // r==3
        auto F = [&](double x,double y){
            vector<double> lam(3);
            lam[0]=x; lam[1]=y; lam[2]=1.0-x-y;
            return evalH1(lam);
        };
        auto inner = [&](double y){
            double l=0.0, rgt = max(0.0, 1.0 - y);
            const double gr = (sqrt(5.0)-1.0)/2.0;
            double x1 = rgt - gr*(rgt - l);
            double x2 = l + gr*(rgt - l);
            double f1 = F(x1,y), f2 = F(x2,y);
            for(int it=0; it<80; ++it){
                if(f1 > f2){
                    l = x1;
                    x1 = x2;
                    f1 = f2;
                    x2 = l + gr*(rgt - l);
                    f2 = F(x2,y);
                } else {
                    rgt = x2;
                    x2 = x1;
                    f2 = f1;
                    x1 = rgt - gr*(rgt - l);
                    f1 = F(x1,y);
                }
            }
            return min(f1,f2);
        };
        double l=0.0, rgt=1.0;
        const double gr = (sqrt(5.0)-1.0)/2.0;
        double y1 = rgt - gr*(rgt - l);
        double y2 = l + gr*(rgt - l);
        double g1 = inner(y1), g2 = inner(y2);
        for(int it=0; it<80; ++it){
            if(g1 > g2){
                l = y1;
                y1 = y2;
                g1 = g2;
                y2 = l + gr*(rgt - l);
                g2 = inner(y2);
            } else {
                rgt = y2;
                y2 = y1;
                g2 = g1;
                y1 = rgt - gr*(rgt - l);
                g1 = inner(y1);
            }
        }
        double ans = min(g1, g2)*100.0;
        cout<<ans<<"\n";
        return 0;
    }
}