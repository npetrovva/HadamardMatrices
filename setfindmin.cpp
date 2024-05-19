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


void WriteToFile(set<vector<uint64_t>> minM, string file_write) { //записываем в файл все Q-эквивив и H-неэквив матрицы для исходной
    ofstream out;
    if (file_write != "\0") {
        out.open(file_write);
    } else {
        out.open("MIN28-1.txt", ios::app);
    }
    for (auto vec: minM) {
        //out << "Q-equiv and H-unequiv matrix is:" << endl;
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

void swap_col0(int c, int s) {
    int tmp,tt0;
    for (int i=0;i<N;i++) {
        if ( H[i] & (((uint64_t)1)<<c) ) { 
            tmp=1; H[i]^=(((uint64_t)1)<<c); 
         }
        else 
            tmp=0;
        if ( H[i] & (((uint64_t)1)<<(s)) ) {
            tt0=1; H[i]^=(((uint64_t)1)<<(s)); 
        }
        else
            tt0=0;
        if (tt0) H[i]^=(((uint64_t)1)<<c);
        if (tmp) H[i]^=(((uint64_t)1)<<(s));
    }
    return;
}

void Column_sort() {
    uint64_t tmp = 1;
    for (int j = 0; j < N - 1; j++) {
        for (int g = j + 1; g < N; g++) {
            bool exch = false;
            int i = 0;
            while (i < N) {
                //int t1 = H[i] & (tmp << (N - 1 - j));
                //int t2 = H[i] & (tmp << (N - 1 - g));
                //cout << t1 << endl;
                if (((H[i] & (tmp << (N - 1 - j))) >> (N - 1 - j)) > (((H[i] & (tmp << (N - 1 - g))) >> (N - 1 - g)))) {
                    swap_col0(N - 1 - j, N - 1 - g);
                    exch = true;
                    break;
                } else if (((H[i] & (tmp << (N - 1 - j))) >> (N - 1 - j)) < (((H[i] & (tmp << (N - 1 - g))) >> (N - 1 - g)))) {
                    break;
                }
                i++;
            }
            if (exch == true) {
                g = j;
            }
            //cout << "j " << j << endl;
            //cout << "g " << g << endl;
        }
    }
    return;
}

void Inversion(uint64_t H_0[N]) {
    uint64_t mask = 1;
    for (int i = 0; i < N; i++) {
        H_0[i] ^= (mask << N) - 1;
    }
    return;
}

void Normalization(uint64_t arr[N]) {
    uint64_t mask = 1;
    for (int i = 0; i < N; i++) { //for 1st column
        //cout << arr[i] << endl;
        //int s = N - 1 - i;
        if (((arr[i] & ((uint64_t)1 << (N - 1))) >> (N - 1)) == 1) {
            arr[i] ^= (((uint64_t)1 << N) - 1);
        }
    }
    
    for (int i = 0; i < N; i++) { //for 1st row
        if (((arr[0] & mask) >> i) == 1) {
            for (int j = 0; j < N; j++) {
                arr[j] ^= ((uint64_t)1) << i;
                //cout << arr[j] <<endl;
            }
        }
        mask <<= 1;
    }
    
    return;
}

void SwapColumns(int i, int j) {
    int tmp1 = 0, tmp2 = 0;
    for (i; i < N; i++) {
        if (H[i] & (((uint64_t)1) << j)) {
            tmp1 = 1;
            H[i] ^= (((uint64_t)1) << j);
        }
        
        if (H[i] & (((uint64_t)1) << (N - 1))) {
            tmp2 = 1;
            H[i] ^= (((uint64_t)1) << (N - 1));
        }

        if (tmp1)
            H[i] ^= (((uint64_t)1) << j);
        if (tmp2)
            H[i] ^= (((uint64_t)1) << (N - 1));
    }

    return;
}

void Core(int r, bool flag) {
    uint64_t tmp = 0;
    if (r == N - 1) { //1
        Column_sort();
        if ((flag == true) or (H[r] < A[r])) //1.2
            A[r] = H[r];
        return; //1.3
    }
    
    uint64_t m  = (((uint64_t)1) << N) - 1; //2
    int k = -1; //3
    int RC[N * 2]; 
    for (int i = r; i < N; i++) { //N or N-1?
        tmp = H[i]; // 4.1
        H[i] = H[r]; //4.1
        H[r] = tmp; //4.1
        Column_sort(); //4.2
        if (H[r] == m) { //4.3
            k++;
            RC[k] = i;
        }
        if (H[r] < m) { //4.4
            k = 0;
            RC[k] = i;
            m = H[r];
        }
        tmp = H[i]; //4.5 re-swap
        H[i] = H[r]; //4.5
        H[r] = tmp; //4.5
    }
    
    if ((flag == 1) or (m < A[r])) { //5
        A[r] = m; //5.1
        tmp = H[r]; // 5.2
        H[r] = H[RC[0]]; //5.2
        H[RC[0]] = tmp; // 5.2
        Column_sort(); //5.3
        Core(r + 1, true); //5.4
        tmp = H[r];//5.5 re-swap
        H[r] = H[RC[0]];//5.5
        H[RC[0]] = tmp; //5.5

        for (int i = 1; i <= k; i++) { //5.6
            tmp = H[r]; //5.6.1
            H[r] = H[RC[i]]; //5.6.1
            H[RC[i]] = tmp; //5.6.1
            Column_sort(); //5.6.2
            Core(r + 1, false); //5.6.3
            tmp = H[r]; //5.6.4 re-swap
            H[r] = H[RC[i]]; //5.6.4
            H[RC[i]] = tmp; //5.6.4
        }
    }

    if ((flag == false) and (m == A[r])) { //6
        for (int i = 0; i <= k; i++) { //6.1
            tmp = H[r]; //6.1.1
            H[r] = H[RC[i]]; //6.1.1
            H[RC[i]] = tmp; //6.1.1
            Column_sort(); //6.1.2
            Core(r + 1, false); //6.1.3
            tmp = H[r]; //6.1.4 re-swap
            H[r] = H[RC[i]]; //6.1.4
            H[RC[i]] = tmp; //6.1.4
        }
    }
    return;
}

void Print_mtrx(uint64_t A[N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int s = N - 1 - j;
            cout << ((A[i] & ((uint64_t)1 << s)) >> s) << ' ';
            if (j == N - 1)
                cout << endl;
        }
    }
    return;
}

void MainAlgorithm() {
    uint64_t tmp = 0;
    memmove(arr, H_0, sizeof(H_0));
    Normalization(arr); //A = Normalization(H_0); //1
    memmove(A, arr, sizeof(arr));
    cout << "Normalized:" << endl;   
    //for (int i = 0; i < N; i++) {
    //    cout << A[i] << endl;
    //}
    //cout << endl;
    Print_mtrx(A);

    memmove(H, H_0, sizeof(H_0)); //2
    
    for (int j = 0; j < N; j++) { //3, swap columns
        //for (int w = 0; w < M; w++) {
          //  swap(H[w][0], H[w][j]); //3.1: swap 1st and j_th columns
	    //}
        SwapColumns(N - 1, N - 1 - j); //3.1
	    for (int i = 0; i < N; i++) { //3.2
            tmp = H[0];
            H[0] = H[i];
            H[i] = tmp; //3.2.1: swap 1st and i_th rows
            memmove(arr, H, sizeof(H));
            Normalization(arr); //H = Normalization(H; //3.2.2
            memmove(H, arr, sizeof(arr));
            Core(1, false); //3.2.3
            tmp = H[0];
            H[0] = H[i];
		    H[i] = tmp; //re-swap 1st and i_th rows
	    }
    }       
    
    for (int i = 0; i < N; i++) 
        H[i] = H_0[i]; //3.3*/
    
    return; //4 return A
}

int main(int argc, char *argv[]) {
    //uint64_t H_0[N];
    set<vector<uint64_t>> minM;
    int flag = 0;
    uint64_t min_matrix[N];
    string line;
    ifstream in;
    //in.open("28-1.txt");
    string filename = argv[1];
    in.open(filename);
    string file_write = "\0";
    if (argc == 3) {
        file_write = argv[2];
    }
    int nus = 0;
    while (getline(in, line)) {
        if (line == "New matrix is:" or line == "Q-equiv matrix is:") {
            for (int y = 0; y < N; y++) {
                H_0[y] = 0;
                cout << H_0[y] << endl;
            }
            nus++;
            cout << "nus " << nus << endl;
            string strcol;
            getline(in, strcol);
            for (int i = 0; i < N; i++) {
                string s;
                getline(in, s);
                cout << "str is: " << s << endl;
                //H_prep[i] = stoull(s);
                //uint64_t t = (uint64_t)1;
                for (int j = N - 1; j >= 0; j--) {
                    if (s[j] == '1') {
                        H_0[i] += pow(2, (N - 1) - j);
                    }
                }
                    
            }
            //flag = 0;

            /*for (int i = 0; i < N; i++) {
                if (!flag) {
                    //H_0[i] = 0;
                    for (int j = 0; j < N + 1; j++) {
                        int a = in.get();
                        if (in.good()) {
                            if (a == 10) {
                                break;
                            }
                            a -= 48;
                            if (a == 1) {
                                H_0[i] |= (uint64_t)1;
                                if (j < N - 1) {
                                    H_0[i] <<= 1;
                                }
                            } else {
                                H_0[i] |= (uint64_t)0;
                                if (j < N - 1) {
                                    H_0[i] <<= 1;
                                }
                            }

                        } else {
                            flag = 1;
                            break;
                        }
                    }   
                }
            }*/
            
            cout << "Matrix H0: " << endl;
            for (int i = 0; i < N; i++) {
                cout << H_0[i] << endl;
            }
            cout << endl;
            /*
            cout << "Input matrix: " << endl;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    int s = N - 1 - j;
                    cout << ((H_0[i] & ((uint64_t)1 << s)) >> s) << ' ';
                    if (j == N - 1)
                        cout << endl;
                }
            }
            cout << endl;*/

            Inversion(H_0);
            
            /*cout << "Inversion Matrix H0: " << endl;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    int s = N - 1 - j;
                    cout << ((H_0[i] & ((uint64_t)1 << s)) >> s) << ' ';
                    if (j == N - 1)
                        cout << endl;
                }
            }
            cout << endl;*/
            
            MainAlgorithm(); //A = MainAlgorithm(H_0)
            
            /*
            cout << "Result: " << endl;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    int s = N - 1 - j;
                    cout << ((A[i] & ((uint64_t)1 << s)) >> s) << ' ';
                    if (j == N - 1)
                        cout << endl;
                }
            }
            cout << endl;*/
            
            vector<uint64_t> stmin;
            for (int i = 0; i < N; i++) {
                stmin.push_back(A[i]);
            }
            minM.insert(stmin);
            //Print_mtrx(A);
            
        }
    }
    in.close();

    WriteToFile(minM, file_write); //запписываем в файл все Q-эквивив и H-неэквив матрицы для исходной

    return 0;
}   
