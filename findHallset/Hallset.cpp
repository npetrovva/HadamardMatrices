#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <stdint.h>
#include <cstring> 
#include <deque>
#include <cmath>
#include <set>
#include <typeinfo>

using namespace std;

#define N 20
uint64_t H_0[N] = {0};
uint64_t qmatrx[N] = {0};
uint64_t arr[N] = {0};
vector<int> strcol;

void WriteToFile(deque<uint64_t> dqFin2, int counter) {
    ofstream out;
    out.open("Rres20-1.txt", ios::app);
    out << "Q-equiv matrix is:" << endl;
    for (auto it = strcol.begin() + counter; it < strcol.begin() + counter + 8; it++) {
        out << *it << ' ';
    }
    out << endl;
    for (auto q : dqFin2) {
        int k = N - 1;
        uint64_t mask = (uint64_t)1 << k;
        while (mask) {
            out << ((q & mask) >> k);
            k--;
            mask >>= 1;
        }
        out << endl;
    }
    out.close();   
    return;
}

void ColumnToString(uint64_t *array, uint64_t *arrbyCol) { 
    int t = N - 1;
    for (int i = 0; i < N; i++) {
        uint64_t tmp = 0, u = 0;
        int k = N - 1;
        for (int j = 0; j < N; j++) {
            u = (array[j] & ((uint64_t)1 << t)) >> t;
            tmp += (u << k);
            k--;
        }
        arrbyCol[i] = tmp;
        t -= 1;
    }
    return;
}

void SwitchHS(uint64_t *arr) { //операция переключения с помощью отрицания F1 и G1
    //int i = 4;
    int j = (N - 4) / 4;
    uint64_t neg = ((uint64_t)1 << (N - 4)) - 1; //00001111111111111111
    cout << "neg = " << neg << endl;
    for (int i = 4; i < 4 + j; i++) {
        arr[i] &= neg;
    }
    return;
}

bool FindG4(deque<uint64_t>& dq) {
    uint64_t value = (((uint64_t)3 << (N - 2))); //1100000...000
    uint64_t value2 = (uint64_t)16 - 1; //1111
    uint64_t negotiation = ((uint64_t)1 << N) - 1; //1111...1111
    uint64_t tmp;
    int i = 4 + 3 * ((N - 4) / 4);
    int j = i;
    for (auto iter = dq.begin() + i; iter < dq.end(); iter++) {
        if ((((*iter) xor value) >> (N - 4)) == (uint64_t)0) {
            tmp = *iter;
            if ((iter - dq.begin()) != i) {
                dq.erase(iter);
                iter--;
                dq.insert(dq.begin() + i, tmp);
            }
            i++;
        } else if ((((*iter) xor value) >> (N - 4)) == value2) {
            tmp = *iter xor negotiation;
            //cout << "HERE2" <<endl;
            if ((iter - dq.begin()) != i) {
                dq.erase(iter);
                iter--;
                dq.insert(dq.begin() + i, tmp);
            } else {
                dq.insert(iter, tmp);
                //cout << "insrt: " << *(iter - 1) << endl;
                //cout << *iter << endl;
                dq.erase(iter + 1);
                //cout << iter - dq.begin() << endl;
            }
            i++;
        }
    }
    cout << "G4 = " << i - j << endl;
    if (i - j == ((N - 4) / 4)) {
        return true;
    } else {
        return false;
    }
}

bool FindG3(deque<uint64_t>& dq) {
    uint64_t value = (((uint64_t)5 << (N - 4))); //0101000...000
    uint64_t value2 = (uint64_t)16 - 1; //1111
    uint64_t negotiation = ((uint64_t)1 << N) - 1; //1111...1111
    uint64_t tmp;
    int i = 4 + 2 * ((N - 4) / 4);
    int j = i;
    for (auto iter = dq.begin() + i; iter < dq.end(); iter++) {
        if (((*iter) xor value) >> (N - 4) == (uint64_t)0) {
            tmp = *iter;
            if ((iter - dq.begin()) != i) {
                dq.erase(iter);
                iter--;
                dq.insert(dq.begin() + i, tmp);
            }
            i++;
        } else if (((*iter) xor value) >> (N - 4) == value2) {
            tmp = *iter xor negotiation;
            if ((iter - dq.begin()) != i) {
                dq.erase(iter);
                iter--;
                dq.insert(dq.begin() + i, tmp);
            } else {
                dq.insert(iter, tmp);
                dq.erase(iter + 1);
            }
            i++;
        }
    }
    cout << "G3 = " << i - j << endl;
    if (i - j == ((N - 4) / 4)) {
        return true;
    } else {
        return false;
    }
}

bool FindG2(deque<uint64_t>& dq) {
    uint64_t value = (((uint64_t)3 << (N - 3))); //0110000...000
    uint64_t value2 = (uint64_t)16 - 1; //1111
    uint64_t negotiation = ((uint64_t)1 << N) - 1; //1111...1111
    uint64_t tmp;
    int i = 4 + (N - 4) / 4;
    int j = i;
    for (auto iter = dq.begin() + i; iter < dq.end(); iter++) {
        if ((((*iter) xor value) >> (N - 4)) == (uint64_t)0) {
            tmp = *iter;
            if ((iter - dq.begin()) != i) {
                dq.erase(iter);
                iter--;
                dq.insert(dq.begin() + i, tmp);
            }
            i++;
        } else if (((*iter) xor value) >> (N - 4) == value2) {
            tmp = *iter xor negotiation;
            if ((iter - dq.begin()) != i) {
                dq.erase(iter);
                iter--;
                dq.insert(dq.begin() + i, tmp);
            } else {
                dq.insert(iter, tmp);
                dq.erase(iter + 1);
            }
            i++;
        }
    }
    cout << "G2 = " << i - j << endl;
    if (i - j == ((N - 4) / 4)) {
        return true;
    } else {
        return false;
    }
}

bool FindG1(deque<uint64_t>& dq) {
    uint64_t value = (((uint64_t)1 << 4) - 1) << (N - 4); //1111000...000
    uint64_t value2 = (uint64_t)16 - 1; //1111
    uint64_t negotiation = ((uint64_t)1 << N) - 1; //1111...1111
    uint64_t tmp;
    int i = 4;
    int j = i;
    for (auto iter = dq.begin() + i; iter < dq.end(); iter++) {
        if (((*iter) & value) == value) {
            tmp = *iter;
            if ((iter - dq.begin()) != i) {
                //cout << (iter - dq.begin()) << endl;
                dq.erase(iter);
                //cout << (iter - dq.begin()) << endl;
                iter--;
                //cout << (iter - dq.begin()) << endl;
                dq.insert(dq.begin() + i, tmp);
                //cout << (iter - dq.begin()) << endl;
            }
            i++;
        } else if ((((*iter) xor value) >> (N - 4)) == value2) {
            tmp = *iter xor negotiation;
            if ((iter - dq.begin()) != i) {
                dq.erase(iter);
                iter--;
                dq.insert(dq.begin() + i, tmp);
            } else {
                //cout << "HERE " <<  iter - dq.begin() << ' ' << endl;
                dq.insert(iter, tmp);
                //cout << (iter - dq.begin()) << endl;
                //cout << *iter << endl;
                //cout << *(iter - 1) << endl;
                dq.erase(iter); // ?????????????? +1
            }
            i++;
        }
    }
    //cout << "G1 = " << i - j << endl;
    if (i - j == ((N - 4) / 4)) {
        return true;
    } else {
        return false;
    }
}

/*
deque<uint64_t> BuildHSMatrix(deque<uint64_t>& HSmatrix, vector<int>& nums, uint64_t *array) {
    for (int i = 0; i < N; i++) {
        auto iter1 { nums.begin() };
        auto iter2 { nums.end() };
        auto res { find(iter1, iter2, i) }; 
        if (res == iter2) { //если не нашли строки в nums
            HSmatrix.push_back(array[i]);
        } else {
            continue;
        }
    }
    return HSmatrix;
}*/
            
/*
set<deque<uint64_t>> FindHallset(deque<uint64_t>& matricesH1, uint64_t *array, deque<vector<int>>& numStr) {
    uint64_t candidate = 0, value = (((uint64_t)1 << 4) - 1) << (N - 4); //1111000...000
    int num = 0;
    deque<uint64_t> HSmatrix;
    set<deque<uint64_t>> result;
    for (auto str: matricesH1) {
        candidate = candidate xor str;
        num++;
        if ((num % 4 == 0) and (candidate == value)) {
            for (int i = num - 4; i < num; i++) {
                HSmatrix.push_back(matricesH1[i]);
            }
            result.insert(BuildHSMatrix(HSmatrix, numStr[num - 4], array));
            candidate = 0;
        } else if (num % 4 == 0) {
            candidate = 0;
            continue;
        }
    }
    return result;
}*/

void FindH4PermCol(uint64_t *array, uint64_t *fixStr, vector<int>& numStr, set<deque<uint64_t>>& res) {
    int idx = 4;
    //дополнили матрицу fixStr
    for (int i = 0; i < N; i++) {
        auto iter1 { numStr.begin() };
        auto iter2 { numStr.end() };
        auto res { find(iter1, iter2, i) }; 
        if (res == iter2) { //если не нашли строки в nums
            fixStr[idx++] = array[i];
            //cout << array[i] << " " ;
        } else {
            continue;
        }
    }
    uint64_t testMatrix[N] = {0};
    ColumnToString(fixStr, testMatrix);
    //cout << endl;
    //перебор столбцов
    //deque<vector<int>> FindH4(uint64_t *array, deque<uint64_t>& newmatrix) {
    //deque<vector<int>> matrixSquare;
    uint64_t m1 = (uint64_t)1 << (N - 1); //1000..
    uint64_t m2 = (uint64_t)1 << (N - 2); //0100..
    uint64_t m3 = (uint64_t)1 << (N - 3); //0010..
    uint64_t m4 = (uint64_t)1 << (N - 4); //0001..
    for (int i = 0; i < N; i++) {
        if ((((testMatrix[i] xor m1) >> (N - 4)) == (uint64_t)0) || (((testMatrix[i] xor m1) >> (N - 4) == (uint64_t)(15)))) {
            for (int j = 0; j < N; j++) {
                if (j != i && ((((testMatrix[j] xor m2) >> (N - 4)) == (uint64_t)0) || (((testMatrix[j] xor m2) >> (N - 4) == (uint64_t)(15))))) {
                    for (int k = 0; k < N; k++) {
                        if (k != i && k != j && ((((testMatrix[k] xor m3) >> (N - 4)) == (uint64_t)0) || (((testMatrix[k] xor m3) >> (N - 4) == (uint64_t)(15))))) {
                            for (int l = 0; l < N; l++) {
                                if (l != k && l != j && l != i && ((((testMatrix[l] xor m4) >> (N - 4)) == (uint64_t)0) || (((testMatrix[l] xor m4) >> (N - 4) == (uint64_t)(15))))) {
                                    uint64_t newmatrix[N] = {0}, findMatrix[N] = {0};
                                    vector<int> cols;
                                    if (((testMatrix[i] xor m1) >> (N - 4)) == (uint64_t)0) {
                                        newmatrix[0] = testMatrix[i];
                                    } else {
                                        newmatrix[0] = testMatrix[i] xor (uint64_t)((1 << N) - 1);
                                    }
                                    if (((testMatrix[j] xor m2) >> (N - 4)) == (uint64_t)0) {
                                        newmatrix[1] = testMatrix[j];
                                    } else {
                                        newmatrix[1] = testMatrix[j] xor (uint64_t)((1 << N) - 1);
                                    }
                                    if (((testMatrix[k] xor m3) >> (N - 4)) == (uint64_t)0) {
                                        newmatrix[2] = testMatrix[k];
                                    } else {
                                        newmatrix[2] = testMatrix[k] xor (uint64_t)((1 << N) - 1);
                                    }
                                    if (((testMatrix[l] xor m4) >> (N - 4)) == (uint64_t)0) {
                                        newmatrix[3] = testMatrix[l];
                                    } else {
                                        newmatrix[3] = testMatrix[l] xor (uint64_t)((1 << N) - 1);
                                    }
                                    //cout << "HERE" << endl;
                                    /*
                                    numStr.push_back(i);
                                    numStr.push_back(j);
                                    numStr.push_back(k);
                                    numStr.push_back(l);*/
                                    for (auto e : numStr) {
                                        strcol.push_back(e);
                                    }
                                    
                                    cols.push_back(i);
                                    cols.push_back(j);
                                    cols.push_back(k);
                                    cols.push_back(l);
                                    for (auto e : cols) {
                                        strcol.push_back(e);
                                    }
                                    /*
                                    for (int x = 0; x < N; x++) {
                                        cout << testMatrix[x] << " ";
                                    }
                                    cout << endl;
                                    */
                                    idx = 4;
                                    //дополнили матрицу newmatrix
                                    for (int i = 0; i < N; i++) {
                                        auto iter1 { cols.begin() };
                                        auto iter2 { cols.end() };
                                        auto res { find(iter1, iter2, i) }; 
                                        if (res == iter2) { //если не нашли строки в nums
                                            newmatrix[idx++] = testMatrix[i];
                                            //cout << "tut: " << i << " " ;
                                        } else {
                                            continue;
                                        }
                                    }
                                    //for (int r = 0; r < N; r++) {
                                    //    cout << newmatrix[r] << " ";
                                    //}
                                    //cout << endl;
                                    uint64_t fMatrix[N] = {0};
                                    ColumnToString(newmatrix, fMatrix); //вернулись к строчному представл.
                                    deque<uint64_t> resMatrix;
                                    for (int u = 0; u < N; u++) {
                                        //cout << fMatrix[u] << " ";
                                        resMatrix.push_back(fMatrix[u]);
                                        //resMatrix.push_back(newmatrix[u]);


                                    }
                                    res.insert(resMatrix);
                                    //пока хз
                                    //matrixSquare.push_back(matrixH4num);
                                    //return matrixH4num;
                                      
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}


set<deque<uint64_t>> FindH4ALL(uint64_t *array) {
    uint64_t matrixSquare[N] = {0};
    uint64_t newar[N] = {0};
    set<deque<uint64_t>> result;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N, i != j; j++) {
            for (int k = 0; k < N, k !=i, k != j; k++) {
                for (int l = 0; l < N, l != i, l != j, l != k; l++) {
                    vector<int> matrixH4num = {16, 3, 1, 15}; //исправить
                    //фиксируем 4 строки
                    matrixSquare[0] = array[i];
                    matrixSquare[1] = array[j];
                    matrixSquare[2] = array[k];
                    matrixSquare[3] = array[l];
//cout << array[i] << " " << array[j] << " " << array[k] << " " << array[l] << endl;
                    //запоминаем номера зафиксированных строк 
                    matrixH4num.push_back(i);
                    matrixH4num.push_back(j);
                    matrixH4num.push_back(k);
                    matrixH4num.push_back(l);
                
                    FindH4PermCol(array, matrixSquare, matrixH4num, result);
                }
            }
        }
    }
    /*vector<int> matrixH4num = {16, 3, 1, 15}; //исправить
    matrixSquare[0] = array[16];
    matrixSquare[1] = array[3];
    matrixSquare[2] = array[1];
    matrixSquare[3] = array[15];
    FindH4PermCol(array, matrixSquare, matrixH4num, result);*/
    return result;
}



void Normalization(uint64_t arr[N]) {
    uint64_t mask = 1;
    for (int i = 0; i < N; i++) { //for 1st column
        if (((arr[i] & ((uint64_t)1 << (N - 1))) >> (N - 1)) == 1) {
            arr[i] ^= (((uint64_t)1 << N) - 1);
        }
    }

    for (int i = 0; i < N; i++) { //for 1st row
        if (((arr[0] & mask) >> i) == 1) {
            for (int j = 0; j < N; j++) {
                arr[j] ^= ((uint64_t)1) << i;
            }
        }
        mask <<= 1;
    }

    return;
}

void PrintM(deque<uint64_t> dq) {
    cout << "matrix3: " << endl;
    //kol++;
    //cout << kol << endl;
    //for (auto d : dq) {
        for (auto str: dq) {
            for (int i = N - 1; i >= 0; i--) {
                cout << ((str & ((uint64_t)1 << i)) >> i) << ' '; 
                if (!i) {
                    cout << endl;
                }
            }
        }
        cout << endl;
    //}
    return;
}

int main() {
    int flag = 0;
    ifstream in;
    in.open("20-1.txt");
    for (int i = 0; i < N; i++) {
        if (!flag) {
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
    }
    in.close();


    cout << "Input matrix: " << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int s = N - 1 - j;
            cout << ((H_0[i] & ((uint64_t)1 << s)) >> s) << ' ';
            if (j == N - 1)
                cout << endl;
        }
    }
    cout << endl;

    memmove(arr, H_0, sizeof(H_0));
    Normalization(arr);
    cout << "Normalized:" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int s = N - 1 - j;
            cout << ((arr[i] & ((uint64_t)1 << s)) >> s) << ' ';
            if (j == N - 1)
                cout << endl;
        }
    }
    cout << endl;
        
    
    //поиск H4 полным перебором
    uint64_t arrbyStr[N] = {0};
    memmove(arrbyStr, arr, sizeof(arr));
    set<deque<uint64_t>> MatrixWithH4 = FindH4ALL(arrbyStr); //все матрицы с H4
    //set<deque<uint64_t>> MatrixWithH4 = MakeByStr(arrbyStr);
    if (MatrixWithH4.size() == 0) {
        cout << "EMPTY" << endl;
    }
    //need
    //PrintM(MatrixWithH4);
    int kol = 0;
    int counter = 0; //ищем  стр и стб
    for (auto dq: MatrixWithH4) {
        if (FindG1(dq)) {
            if (FindG2(dq)) {
                if (FindG3(dq)) {
                    if (FindG4(dq)) {
                       

                         //need
                        uint64_t dqCol[N] = {0};
                        uint64_t tmparr[N] = {0};
                        int i = 0;
                        for (auto s: dq) {
                            tmparr[i++] = s;
                        }
                        ColumnToString(tmparr, dqCol);
                        deque<uint64_t> dqFin;
                        //need
                        for (i = 0; i < N; i++) {
                            dqFin.push_back(dqCol[i]);
                        }
                        //need    
                        if (FindG1(dqFin)) {
                            if (FindG2(dqFin)) {
                                if (FindG3(dqFin)) {
                                    if (FindG4(dqFin)) {
                                        i = 4 + (N - 4) / 4;
                                        int e = 0;
                                        deque<uint64_t> dqFin1;
                                        uint64_t neg = ((uint64_t)1 << N) - 1;
                                        for (auto d: dqFin) {
                                            if (e >= i) {
                                                dqFin1.push_back(d xor neg);
                                            } else {
                                                dqFin1.push_back(d);
                                            }
                                            e++;
                                        }

                                        cout << endl;
                                        
                                        i = 0;
                                        for (auto s: dqFin1) {
                                            tmparr[i++] = s;
                                        }
                                        SwitchHS(tmparr); //F1
                                        ColumnToString(tmparr, dqCol);
                                        SwitchHS(dqCol); //G1
                                        deque<uint64_t> dqFin2;
                                        for (i = 0; i < N; i++) {
                                            dqFin2.push_back(dqCol[i]);
                                        }
                                            
                                        cout << "HSmatrix3: " << endl;
                                        kol++;
                                        cout << kol << endl;
                                        for (auto str: dqFin2) {
                                            for (int i = N - 1; i >= 0; i--) {
                                                cout << ((str & ((uint64_t)1 << i)) >> i) << ' '; 
                                                if (!i) {
                                                    cout << endl;
                                                }
                                            }
                                        }
                                        cout << endl;
                                        
                                        WriteToFile(dqFin2, counter);
                                        counter += 8;
                                    } else {
                                        counter += 8;
                                        continue;
                                    }
                                } else {
                                    counter += 8;
                                    continue;
                                }
                            } else {
                                counter += 8;
                                continue;
                            }
                        } else {
                            counter += 8;
                            continue;
                        }
                    } else {
                        counter += 8;
                        continue;
                    }
                } else {
                    counter += 8;
                    continue;
                }
            } else {
                counter += 8;
                continue;
            } 
        } else {
            counter += 8;
            continue;
        }
                        
    }

    return 0;
}
