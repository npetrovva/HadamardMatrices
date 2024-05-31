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

#define N 24
uint64_t H_0[N] = {0};
uint64_t qmatrx[N] = {0};
uint64_t arr[N] = {0};
vector<int> goodcols;

void PrintMatrix(uint64_t *a) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int s = N - 1 - j;
            cout << ((a[i] & ((uint64_t)1 << s)) >> s) << ' ';
            if (j == N - 1)
                cout << endl;
        }
    }
    cout << endl;
    return;
}

void WriteToFile(set<deque<uint64_t>> res, int counter, string file_write) {
    ofstream out;
    if (file_write != "\0") {
        out.open(file_write, ios::app);
    } else {
        out.open("minRres24-1.txt", ios::app);
    }
    for (auto dq: res) {
        out << "Q-equiv matrix is:" << endl;
        for (auto it = goodcols.begin() + counter; it < goodcols.begin() + counter + 4; it++) {
            out << *it << ' ';
        }
        
        out << endl;
        for (auto q : dq) {
            int k = N - 1;
            uint64_t mask = (uint64_t)1 << k;
            while (mask) {
                out << ((q & mask) >> k);
                k--;
                mask >>= 1;
            }
            out << endl;
        }
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


bool ChecktheStruct(uint64_t *matrix, int k, uint64_t *matrix_sort) {
    //matrix with first 4 str: 0,t,j,k
    uint64_t matrix_col[N];
    //uint64_t matrix_sort[N];
    int u, v;
    int l = 0, m = 0;
    memset(matrix_col, 0, N * sizeof(uint64_t));
    memset(matrix_sort, 0, N * sizeof(uint64_t));
    memmove(matrix_col, matrix, N * sizeof(uint64_t));
    uint64_t mask = (uint64_t)1 << (N - 2); //010...0
    if (k == 1) {
        u = N - 1, v = 0;
        for (int i = 0; i < N; i++) {
            if (((matrix_col[i] & mask) >> (N - 2)) == 1) {
                matrix_sort[u--] = matrix_col[i];
                l++;
            } else {
                matrix_sort[v++] = matrix_col[i];
                m++;
            }
        }

    
        if (l != m)
            return false;
     
    } else if (k == 2) {
        mask >>= 1;
        u = N/2 - 1, v = 0;
        for (int i = 0; i < N/2; i++) {
            if (((matrix_col[i] & mask) >> (N - 3)) == 1) {
                matrix_sort[u--] = matrix_col[i];
                l++;
            } else {
                matrix_sort[v++] = matrix_col[i];
                m++;
            }
        
        }
        if (l != m && l != 0 && m != 0)
            return false;
        l = 0, m = 0;
      
        u = N - 1, v = N/2;
        for (int i = N/2; i < N; i++) {
            if (((matrix_col[i] & mask) >> (N - 3)) == 1) {
                matrix_sort[u--] = matrix_col[i];
                l++;
            } else {
                matrix_sort[v++] = matrix_col[i];
                m++;
            }

        }
        if (l != m && l != 0 && m != 0) 
            return false;
      
    
    } else if (k == 3) {
        mask >>= 1;
        u = N/4 - 1, v = 0;
        for (int i = 0; i < N/4; i++) {
            if (((matrix_col[i] & mask) >> (N - 4)) == 1) {
                matrix_sort[u--] = matrix_col[i];
                l++;
            } else {
                matrix_sort[v++] = matrix_col[i];
                m++;
            }
        
        }
        /*if (v - u != 1) {
            cout << "3. " << "v = " << v << " " << "u = " << u << endl;
            return false;
        }*/
        if (l > 0 && m > 0)
            return false;
        l = 0, m = 0;
        u = N/2 - 1, v = N/4;
        for (int i = N/2; i < N; i++) {
            if (((matrix_col[i] & mask) >> (N - 4)) == 1) {
                matrix_sort[u--] = matrix_col[i];
                l++;
            } else {
                matrix_sort[v++] = matrix_col[i];
                m++;
            }

        }
      
        if (l > 0 && m > 0)
            return false;
        l = 0, m = 0;
        u = 3 * N/4 - 1, v = N/2;
        for (int i = 0; i < N/4; i++) {
            if (((matrix_col[i] & mask) >> (N - 4)) == 1) {
                matrix_sort[u--] = matrix_col[i];
                l++;
            } else {
                matrix_sort[v++] = matrix_col[i];
                m++;
            }
        
        }
        /*if (v - u != 1) {
            return false;
        }*/
        if (l > 0 && m > 0)
            return false;
        l = 0, m = 0;

        u = N - 1, v = 3 * N/4;
        for (int i = N/2; i < N; i++) {
            if (((matrix_col[i] & mask) >> (N - 4)) == 1) {
                matrix_sort[u--] = matrix_col[i];
                l++;
            } else {
                matrix_sort[v++] = matrix_col[i];
                m++;
            }

        }
        if (l > 0 && m > 0)
            return false;
       
    }
    return true;

}

set<deque<uint64_t>>FindQuadr(uint64_t *matrixWithFixStr, int i) {
    set<deque<uint64_t>> res;
    uint64_t matrixSquare[N] = {0};
     
    cout << "Before norm:" << endl;
    //PrintMatrix(matrixWithFixStr);
    
    //cout << "=======================" << endl;
    uint64_t mask = matrixWithFixStr[0];
    uint64_t mask2 = ((uint64_t)1 << N) - 1; //1111...111
    for (int t = 0; t < N; t++) { //for 1st row
        matrixWithFixStr[t] ^= mask;
        matrixWithFixStr[t] ^= mask2;

    }
    /*
    for (int t = 0; t < N; t++) { //for 1st row
        cout << matrixWithFixStr[t] << endl;
    }*/
    //PrintMatrix(matrixWithFixStr);
    //cout << "=======================" << endl;
    /*
    if (i > 0) {
    //если i-я строка ненулевая, производим нормализацию по этой строке путем перестановки столбцов (только строки, без столбца) 
        uint64_t mask = 1;
        for (int t = 0; t < N; t++) { //for 1st row
            if (((matrixWithFixStr[0] & mask) >> t) == 0) {
                for (int j = 0; j < N; j++) {
                    matrixWithFixStr[j] ^= ((uint64_t)1) << t;
                }
            }
            mask <<= 1;
        }
    }
    */
    //ищем 3 строки среди оставшихся N-1
    for (int t = 1; t < N; t++) {
        for (int j = 1; j < N; j++) {
            for (int k = 1; k < N; k++) {
                if (t != j && t != k && j != k) {
                    //проверяем, что выполняется три-нормализация для этих строк
                    if ((matrixWithFixStr[t] xor matrixWithFixStr[j] xor matrixWithFixStr[k]) == mask2) {
                        //cout << "HERE" << endl;
                        //check the structure of quadruple
                        //проверяем, что полученная четверка строк соответствуем необходимой структуре, если нет - ищем новую четверку в цикле
                        
                        /*uint64_t m1 = matrixWithFixStr[t], m2 = matrixWithFixStr[j], m3 = matrixWithFixStr[k]; 
                        if ((ChecktheStruct(m1) == false) or (ChecktheStruct(m2) == false) or (ChecktheStruct(m3) == false)) {
                            continue;
                        }*/
                        uint64_t tmp_m[N];
                        memset(tmp_m, 0, N * sizeof(uint64_t));
                        tmp_m[0] = matrixWithFixStr[0];
                        tmp_m[1] = matrixWithFixStr[t];
                        tmp_m[2] = matrixWithFixStr[j];
                        tmp_m[3] = matrixWithFixStr[k];
                        int a = 4;
                        for (int l = 1; l < N; l++) {
                            if (l != j && l != t && l != k) {
                                tmp_m[a++] = matrixWithFixStr[l];
                            }
                        }
                        
                        //PrintMatrix(tmp_m);
                        uint64_t colstmp[N];
                        memset(colstmp, 0, N * sizeof(uint64_t));
                        ColumnToString(tmp_m, colstmp);
                        uint64_t new_m[N];
                        memset(new_m, 0, N * sizeof(uint64_t));
                        bool fquadr = true;
                        for (int i = 1; i < 3; i++) {
                            if (ChecktheStruct(colstmp, i, new_m) == false) {
                                fquadr = false;
                                break;
                            }
                            memmove(colstmp, new_m, N * sizeof(uint64_t));
                            cout <<"+++++++++++++++++++++++++++++++++++" << endl;
                            PrintMatrix(colstmp);
                        }
                        if (fquadr == false) {
                            cout << "HEREFALSE" << endl;
                            continue;
                        }
                        ColumnToString(colstmp, tmp_m);
                        cout <<"After finfing quadr : " << endl;
                        PrintMatrix(tmp_m);
                        cout << "tjk: " << t << j << k << endl;
                        vector<int> numStr;
                        //перезаписываем матрицу из массива в очередь
                        deque<uint64_t> matrix;
                        for (int u = 0; u < N; u++) {
                            matrix.push_back(tmp_m[u]);
                        }
                        /*
                        matrix.push_back(matrixWithFixStr[0]);
                        matrix.push_back(matrixWithFixStr[t]);
                        matrix.push_back(matrixWithFixStr[j]);
                        matrix.push_back(matrixWithFixStr[k]);
                        for (int u = 1; u < N; u++) {
                            if (u != t && u != j && u != k) {
                                matrix.push_back(matrixWithFixStr[u]);
                            }
                        }*/
         		 //фиксируем 3 строки
                        /*
                        matrixSquare[0] = matrixWithFixStr[0];
                        matrixSquare[1] = matrixWithFixStr[t];
                        matrixSquare[2] = matrixWithFixStr[j];
                        matrixSquare[3] = matrixWithFixStr[k];
                        
                        numStr.push_back(i);
                        numStr.push_back(t);
                        numStr.push_back(j);
                        numStr.push_back(k);
                        */
                        //здесь я хотела зафиксировать номера строк, но они уплыли из-за перемещения i-й строки на 1 ппозицию и перезаписи остальных
                        goodcols.push_back(i);
                        goodcols.push_back(t);
                        goodcols.push_back(j);
                        goodcols.push_back(k);
                        //добавляю матрицу с замкнутой четверкой в сет
                        res.insert(matrix);

                    }
                }
            }
        }
    }
    return res;
}
        

set<deque<uint64_t>> FixFirstStr(uint64_t *array) {
    set<deque<uint64_t>> result;
    uint64_t matrixWithFixStr[N];
    memset(matrixWithFixStr, 0, N * sizeof(uint64_t));
    int i;
    //каждую строку ставим на 1 место для последующей нормализации
    for (int i = 0; i < N; i++) {
        uint64_t first = array[i];
        if (first != array[0]) {
            uint64_t tmp = array[0];
            matrixWithFixStr[0] = first;
            //дополняем матрицу, когда поставили i-ю строку на 1 место
            for (int j = 1; j < N; j++) { 
                if (j > i) {
                    matrixWithFixStr[j] = array[j]; 
                } else {
                    matrixWithFixStr[j] = array[j - 1];
                }
            }
        } else {
            memmove(matrixWithFixStr, array, sizeof(array));
        }
        /*
        for (int i = 0; i < N; i++) {
            cout << matrixWithFixStr[i] << endl;
        }
        cout << "--------------";
        */
        //передаем далее на обработку матрицу с i-й строкой на 1 месте
        set<deque<uint64_t>> resforOneStr = FindQuadr(matrixWithFixStr, i);
        //полученный сет матриц для i-й строки на 1 месте добавляю в общий сет матриц с замкнутыми четвеками для исходной матрицы
        for (auto item: resforOneStr) {
            result.insert(item);
        }
    
    }
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

int main(int argc, char *argv[]) {
    int flag = 0;
    ifstream in;
    //in.open("min24-1.txt");
    string filename = argv[1];
    in.open(filename);
    string file_write = "\0";
    if (argc == 3) {
        file_write = argv[2];
    }
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
            //out << ((arr[i] & ((uint64_t)1 << s)) >> s) << ' ';
            if (j == N - 1) {
                cout << endl;
                //out << endl;
            }
        }
    }
    cout << endl;

    uint64_t arrbyStr[N] = {0};
    //uint64_t arrbyCol[N] = {0};
    memmove(arrbyStr, arr, sizeof(arrbyStr));
    //ColumnToString(arr, arrbyCol);
    set<deque<uint64_t>> MatrixWithQuadr;
    set<deque<uint64_t>> MatrixWithQuadrByStrs = FixFirstStr(arrbyStr);
    if (MatrixWithQuadrByStrs.size() == 0) {
        cout << "There is no closed quadruple strs" << endl;
    } else {
        cout << "Strs is ok" << endl;
        cout << "size: " << MatrixWithQuadrByStrs.size() << endl;
        for (auto t: MatrixWithQuadrByStrs) {
            MatrixWithQuadr.insert(t);
        }
    }
    /*
    set<deque<uint64_t>> MatrixWithQuadrByCols = FixString(arrbyCol);
    if (MatrixWithQuadrByCols.size() == 0) {
        cout << "There is no closed quadruple cols" << endl;
    } else {
        cout << "Cols is ok" << endl;
        for (auto t: MatrixWithQuadrByCols) {
            MatrixWithQuadr.insert(t);
        }
    }
    */

    //построение всех матриц, q-эквивалентных заданной,
    //методом переключения замкнутой четверки
    uint64_t mask = (((uint64_t)1 << (N / 4)) - 1) << (N - (N / 4)); //1111000000
    set<deque<uint64_t>> result;
    int count = 0;
    for (auto firstQuadr: MatrixWithQuadr) { //матрица в сете
        //deque<uint64_t> qmatrx;
        int i = 0;
        for (auto iter = firstQuadr.begin(); iter < firstQuadr.begin() + 4; iter++) {
            auto tmp = *iter;
            i = iter - firstQuadr.begin();
            tmp = tmp xor mask;
            firstQuadr.erase(iter);
            firstQuadr.insert(firstQuadr.begin() + i, tmp);
            iter = firstQuadr.begin() + i;
            i++;
        }
        result.insert(firstQuadr);
        count += 4;
        /*
        cout << endl;
        for (int i = 0; i < N; i++) {
            bool used = false;
            for (int s: firstQuadr) {
                if (i == s) {
                    used = true;
                    break;
                }
            }
            if (!used) {
                qmatrx.push_back(arr[i]);
            }
        }

        for (int i = 0; i < qmatrx.size(); i++) {
            for (int j = 0; j < N; j++) {
                int s = N - 1 - j;
                cout << ((qmatrx[i] & ((uint64_t)1 << s)) >> s) << ' ';
                if (j == N - 1)
                    cout << endl;
            }
        }
        cout << endl;
        result.insert(qmatrx);*/
    }
    WriteToFile(result, count, file_write); //запись полученных матриц в файл
    return 0;

}
