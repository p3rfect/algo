#include <bits/stdc++.h>

#pragma GCC optimize("O3")

#define ll long long
#define INF (ll)1e17

using namespace std;

//Any K <= 23 is good with M = 998244353 or M = 9223372036737335297LL (if input numbers are of long long type, then use __int128)
template<int M, int K, int G> struct FFT{
    int n, A[1 << K], B[1 << K];

    int mod(int a){
        return (a % M + M) % M;
    }

    int bin_pow(int a, int b){
        int res = 1;
        while (b > 0){
            if (b & 1) res = 1LL * res * a % M;
            a = 1LL * a * a % M;
            b >>= 1;
        }
        return res;
    }

    void fft(bool inv){
        fill(B, B + n, 0);
        for (int b = n / 2; b > 0; b /= 2, swap(A, B)){
            int w = bin_pow(G, (M - 1) / n * b), m = 1;
            for (int i = 0; i < n; i += b * 2, m = 1LL * m * w % M){
                for (int j = 0; j < b; j++){
                    int u = A[i + j], v = 1LL * A[i + j + b] * m % M;
                    B[i / 2 + j] = mod(u + v);
                    B[i / 2 + j + n / 2] = mod(u - v);
                }
            }
        }
        if (inv){
            reverse(A + 1, A + n);
            int z = bin_pow(n, M - 2);
            for (int i = 0; i < n; i++){
                A[i] = 1LL * A[i] * z % M;
            }
        }
    }

    int temp[1 << K];

    void init(vector < int > &x){
        for (int i = 0; i < n; i++) A[i] = x[i];
    }

    void multiply(vector < int > &x, vector < int > &y){
        int _n = max(x.size(), y.size()) * 2;
        int val = 1;
        while (val < _n) val *= 2;
        n = val;
        while (x.size() < n) x.push_back(0);
        while (y.size() < n) y.push_back(0);
        init(x);
        fft(0);
        swap(temp, A);
        init(y);
        fft(0);
        for (int i = 0; i < n; i++){
            A[i] = 1LL * A[i] * temp[i] % M;
        }
        fft(1);
    }
};

FFT<998244353,15,31> muls;

int n, i;

vector < int > a, b;

void solve(){
    cin >> n;
    a.resize(n);
    b.resize(n);
    for (i = 0; i < n; i++) cin >> a[i];
    for (i = 0; i < n; i++) cin >> b[i];
    muls.multiply(a, b);

    for (int i = 0; i < n * 2 - 1; i++){
        cout << muls.A[i] << " ";
    }
    cout << '\n';
}

int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    int q;
    cin >> q;

    while (q--){
        solve();
    }
    return 0;
}
