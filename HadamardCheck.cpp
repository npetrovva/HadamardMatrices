#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <cmath>
#include <set>

using namespace std;

#define N 28

uint64_t A[N] = {0};
uint64_t H[N] = {0};
uint64_t H_0[N] = {0};
uint64_t H_prep[N] = {0};
uint64_t arr[N] = {0};
int notH = 0; //сколько неадамаровых
int isH = 0;

void WriteToFile(set<vector<uint64_t>> minM, vector<string> s) { //записываем в файл все Q-эквивив и H-неэквив матрицы для исходной
    ofstream out;
    out.open("NOTRres16-1.txt", ios::app);
    auto iter = s.begin();
    for (auto vec: minM) {
        out << "Q-equiv and H-unequiv matrix is:" << endl;
        //out << "Q-equiv matrix is:" << endl;
        //out << *iter << endl;
        iter++;
        for (auto s : vec) {
            int k = N - 1;
            uint64_t mask = (uint64_t)1 << k;
            while (mask) {
                out << ((s & mask) >> k);
                k--;
                mask >>= 1;
            }
            out << endl;
        }
    }
    out.close();   
    return;
}

bool HadamardCheck(vector<uint64_t>& m) {
    for (auto it1 = m.begin(); it1 < m.end() - 1; it1++) {
        for (auto it2 = it1 + 1; it1 != it2, it2 < m.end(); it2++) {
            uint64_t st = *it1 xor *it2;
            //cout << "st" << st << endl;
            uint64_t mask = (uint64_t) 1;
            int ones = 0;
            int zeroes = 0;
            int n = 28;
            while (n) {
                if (st & mask == 1) {
                    ones++;
                    //cout << "ones " << ones << endl;
                } else {
                    zeroes++;
                    //cout << "zeroes " << zeroes << endl;
                }
                st >>= 1;
                n--;
            }
            //cout << "ones = " << ones << ' ' << "zeroes = " << zeroes << endl;
            if (zeroes != ones) {
                cout << "ones = " << ones << ' ' << "zeroes = " << zeroes << endl;
                cout << "it1 = " << it1 - m.begin() << endl;
                cout << "it2 = " << it2 - m.begin() << endl;
                cout << "*it1 = " << *it1 << endl;
                cout << "*it2 = " << *it2 << endl;
                notH++;
                return false;
            }
        }
    }
    isH++;
    return true;
}

int main() {
    //uint64_t H_0[N];
    set<vector<uint64_t>> isHadamard;
    int flag = 0;
    //uint64_t matrix[N];
    string line;
    ifstream in;
    in.open("minHad.28.34");
    //in.open("Rres16-1.txt");
    int nus = 0;
    vector<string> s;
    while (getline(in, line)) {
        if (line == "Q-equiv matrix is:") {
            //memset(H_0, 0, sizeof(int) * N);
            for (int y = 0; y < N; y++) {
                H_0[y] = 0;
                //cout << H_0[y] << endl;
            }
            nus++;
            cout << endl;
            cout << "nus " << nus << endl;
            string strcol;
            getline(in, strcol);
            for (int i = 0; i < N; i++) {
                string s;
                getline(in, s);
                //cout << "str is: " << s << endl;
                //H_prep[i] = stoull(s);
                //uint64_t t = (uint64_t)1;
                for (int j = N - 1; j >= 0; j--) {
                    if (s[j] == '1') {
                        H_0[i] += pow(2, (N - 1) - j);
                    }
                }

            }
            /*cout << "Matrix H0: " << endl;
            for (int i = 0; i < N; i++) {
                cout << H_0[i] << endl;
            }
            cout << endl;*/

            vector<uint64_t> matrix;
            for (int i = 0; i < N; i++) {
                matrix.push_back(H_0[i]);
            }
            
            if (HadamardCheck(matrix) == false) {
                s.push_back(strcol);
                isHadamard.insert(matrix);
            }
            //Print_mtrx(A);

        }
    }
    in.close();

    WriteToFile(isHadamard, s);
    
    cout << "NOTH = " << notH << endl;
    cout << "ISH = " << isH << endl;
    return 0;
}
