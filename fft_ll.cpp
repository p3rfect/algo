#include <bits/stdc++.h>

#define ll long long

using namespace std;

template<__int128 M, __int128 K, __int128 G> struct FFT{ //any K <= 23 is good with M = 998244353
    __int128 n, A[1 << K], B[1 << K];

    __int128 mod(__int128 a){
        return (a % M + M) % M;
    }

    __int128 bin_pow(__int128 a, __int128 b){
        __int128 res = 1;
        while (b > 0){
            if (b & 1) res = res * a % M;
            a = a * a % M;
            b >>= 1;
        }
        return res;
    }

    void fft(bool inv){
        fill(B, B + n, 0);
        for (int b = n / 2; b > 0; b /= 2, swap(A, B)){
            __int128 w = bin_pow(G, (M - 1) / n * b), m = 1;
            for (int i = 0; i < n; i += b * 2, m = m * w % M){
                for (int j = 0; j < b; j++){
                    __int128 u = A[i + j], v = A[i + j + b] * m % M;
                    B[i / 2 + j] = mod(u + v);
                    B[i / 2 + j + n / 2] = mod(u - v);
                }
            }
        }
        if (inv){
            reverse(A + 1, A + n);
            __int128 z = bin_pow(n, M - 2);
            for (int i = 0; i < n; i++){
                A[i] = A[i] * z % M;
            }
        }
    }

    __int128 temp[1 << K];

    void init(vector < ll > &x){
        for (int i = 0; i < n; i++) A[i] = x[i];
    }

    void multiply(vector < ll > &x, vector < ll > &y){
        __int128 _n = max(x.size(), y.size()) * 2;
        __int128 val = 1;
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

FFT<9223372036737335297LL,15,31> muls;

ll n, i;

vector < ll > a, b;

void solve(){
    cin >> n;
    a.resize(n + 1);
    b.resize(n + 1);
    for (i = 0; i < n; i++) cin >> a[i];
    for (i = 0; i < n; i++) cin >> b[i];
    muls.multiply(a, b);

    for (ll i = 0; i < n * 2; i++){
        cout << (ll)muls.A[i] << " ";
    }
    cout << '\n';
}

int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    ll q;
    cin >> q;

    while (q--){
        solve();
    }
    return 0;
}
