//
//  main.cpp
//  KnightTour Sequential
//
//  Created by Serguei Diaz on 21.12.2024.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <chrono>

#include "baseSolutions.hpp"

using namespace std;

enum Corner {
    topLeft,
    topRight,
    bottomLeft,
    bottomRight,
    none
};

enum TourType {
    open,
    stretched
};

enum JoinType {
    horizontal,
    horizontalDoubleLoop,
    horizontalStretchedDoubleLoop,
    vertical
};

enum stretchedType {
    horizontalStretched,
    verticalStretched,
    anyStretched,
    noStretched
};

vector<vector<int>> optimalOpenKnight(int N, int M, stretchedType stretchedType);
vector<vector<int>> secondCaseSolution(int n, int m);
vector<vector<int>> join(vector<vector<int>> a, vector<vector<int>> b, int n1, int m1, int n2, int m2, JoinType type );
void reverseMatrixValues(vector<vector<int>>& matrix, int maxValue);
void addToMatrixConditional(vector<vector<int>>& matrix, int value, int startFrom);
void invertMatrixHorizontally(vector<vector<int>>& matrix);

void printBoard(const vector<vector<int>>& board) {
    for (const auto& row : board) {
        cout << "{ ";
        for (int cell : row) {
            cout << cell << ", ";
        }
        cout << " },";
        cout << endl;
    }
}

void reverseMatrixValues(vector<vector<int>>& matrix, int maxValue) {
    for (size_t row = 0; row < matrix.size(); ++row) {
        for (size_t column = 0; column < matrix[row].size(); ++column) {
            matrix[row][column] = maxValue - matrix[row][column];
        }
    }
}

void minusMatrixValues(vector<vector<int>>& matrix) {
    for (size_t row = 0; row < matrix.size(); ++row) {
        for (size_t column = 0; column < matrix[row].size(); ++column) {
            matrix[row][column] = - matrix[row][column] - 2;
        }
    }
}

void addToMatrixConditional(vector<vector<int>>& matrix, int value, int startFrom) {
    for (size_t row = 0; row < matrix.size(); ++row) {
        for (size_t column = 0; column < matrix[row].size(); ++column) {
            if (matrix[row][column] >= startFrom) {
                matrix[row][column] = matrix[row][column] - startFrom + value;
            }
        }
    }
}

int addToMatrixDoubleConditional(vector<vector<int>>& matrix, int value, int startFrom, int until, int mult, bool reverse) {
    int minMax = -1;
    for (size_t row = 0; row < matrix.size(); ++row) {
        for (size_t column = 0; column < matrix[row].size(); ++column) {
            if (matrix[row][column] >= startFrom && matrix[row][column] <= until) {
                if (reverse) {
                    matrix[row][column] = mult * (matrix[row][column] - until) + value;
                }
                else {
                    matrix[row][column] = mult * (matrix[row][column] - startFrom) + value;
                }
                
                if (minMax == -1) {
                    minMax = matrix[row][column];
                }
                if (mult > 0) {
                    if (minMax < matrix[row][column]) {
                        minMax = matrix[row][column];
                    }
                }
                else {
                    if (minMax > matrix[row][column]) {
                        minMax = matrix[row][column];
                    }
                }
                
            }
        }
    }
    
    return minMax;
}

void invertMatrixHorizontally(vector<vector<int>>& matrix) {
    for (auto& row : matrix) {
        reverse(row.begin(), row.end());
    }
}

vector<vector<int>> joinMatricesHorizontally(const vector<vector<int>> matrix1, const vector<vector<int>> matrix2) {
    if (matrix1.size() != matrix2.size()) {
        throw std::invalid_argument("Matrices must have the same number of rows to join horizontally.");
    }
    
    vector<vector<int>> result = matrix1;
    for (size_t i = 0; i < matrix1.size(); ++i) {
        result[i].insert(result[i].end(), matrix2[i].begin(), matrix2[i].end());
    }
    
    return result;
}

vector<vector<int>> joinMatricesVertically(const vector<vector<int>>& matrix1, const vector<vector<int>>& matrix2) {
    if (!matrix1.empty() && !matrix2.empty() && matrix1[0].size() != matrix2[0].size()) {
        throw std::invalid_argument("Matrices must have the same number of columns to join vertically.");
    }

    vector<vector<int>> result = matrix1;
    result.insert(result.end(), matrix2.begin(), matrix2.end());

    return result;
}

vector<vector<int>> join(vector<vector<int>> a, vector<vector<int>> b, int n1, int m1, int n2, int m2, JoinType type ) {
    if (a.empty()) {
        return b;
    }
    else if (b.empty()) {
        return a;
    }
    
    int temp;
    int valTemp;
    int valsecondTemp;
    
    switch (type) {
        case horizontal:
            if (a[2][m1-1] > a[0][m1-2]) {
                reverseMatrixValues(b, (n2 * m2 - 1));
                addToMatrixConditional(b, a[0][m1-2] + 1, 0);
                addToMatrixConditional(a, b[0][0] + 1, a[2][m1-1]);
                return joinMatricesHorizontally(a, b);
            }
            else {
                addToMatrixConditional(b, a[2][m1-1] + 1, 0);
                addToMatrixConditional(a, b[1][0] + 1, a[0][m1-2]);
                return joinMatricesHorizontally(a, b);
            }
            
        case vertical:
            if (a[n1-1][2] > a[n1-2][0]) {
                reverseMatrixValues(b, (n2 * m2 - 1));
                addToMatrixConditional(b, a[n1-2][0] + 1, 0);
                addToMatrixConditional(a, b[0][0] + 1, a[n1-1][2]);
                return joinMatricesVertically(a, b);
            }
            else {
                addToMatrixConditional(b, a[n1-1][2] + 1, 0);
                addToMatrixConditional(a, b[0][1] + 1, a[n1-2][0]);
                addToMatrixConditional(a, b[0][1] + 1, a[n1-2][0]);
                return joinMatricesVertically(a, b);
            }
            
        case horizontalStretchedDoubleLoop:
            invertMatrixHorizontally(a);
            
            addToMatrixConditional(a, (m2 * n2 / 2), 0);
            
            valTemp = b[1][0] - 1;
            temp = addToMatrixDoubleConditional(b, a[0][m1-2] + 1, b[1][0], (m2 * n2 / 2) - 1, 1, false);
            temp = addToMatrixDoubleConditional(b, temp + 1, 0, valTemp, 1, false);
            
            valsecondTemp = b[2][0] - 1 + (m2 * n2 / 2) + 2;
            temp = addToMatrixDoubleConditional(b, 0, -(m2 * n2 / 2) - 1, b[2][0] - 1, -1, true);
            temp = addToMatrixDoubleConditional(b, valsecondTemp, b[2][0], -2, -1, true);
            
            return joinMatricesHorizontally(a, b);
            
        case horizontalDoubleLoop:
            minusMatrixValues(b);
            
            valTemp = a[n1-2][m1-1] - 1;
            valsecondTemp = a[n1-2][m1-1] + 2 + (b[n2-1][1] - 1) - 1;

            temp = addToMatrixDoubleConditional(a, valsecondTemp, -(m1 * n1 / 2) - 1, valTemp, -1, false);
            
            addToMatrixDoubleConditional(a, b[n2-1][1] - 1, a[n1-2][m1-1], -1, -1, false);
            
            temp = addToMatrixDoubleConditional(b, temp - 1, -(m2 * n2 / 2) - 1, b[1][0], 1, true);
            
            valsecondTemp = a[1][m1-1] + 1 + b[0][1] + 1;
            valTemp = b[0][1] + 1;
            
            temp = addToMatrixDoubleConditional(b, valsecondTemp, valTemp, (m2 * n2 / 2) - 1, -1, true);
            
            
            valsecondTemp = a[1][m1-1] + 1 + b[0][1] + 1;

            temp = addToMatrixDoubleConditional(b, a[1][m1-1] + 1, 0, b[0][1], -1, true);
            
            temp = addToMatrixDoubleConditional(a, b[2][0] + 1, a[n1-1][m1-2], (m2 * n2 / 2), 1, false);
            
            return joinMatricesHorizontally(a, b);
    }
}

vector<vector<int>> firstCaseSolution(int n, int m, stretchedType stretchedType) {
    if (stretchedType == horizontalStretched) {
        return getStretchedSolutionH(n, m);
    }
    else if (stretchedType == verticalStretched) {
        return getStretchedSolutionV(n, m);
    }
    else {
        if ( !getOpenSolution(n, m).empty()){
            return getOpenSolution(n, m);
        }
        else {
            vector<vector<int>> response = getStretchedSolutionH(n, m);
            if (response.empty()) {
                response = getStretchedSolutionV(n, m);
            }
            return response;
        }
    }
}

vector<vector<int>> secondCaseSolution(int n, int m) {
    if (n != 3 || m <= 10) {
        return vector<vector<int>>();
    }
    
    int k = ((m - 7) % 4) + 7;
    
    vector<vector<int>> left = getOpenSolution(n, k);
    
    vector<vector<int>> temp0;
    
    vector<vector<int>> temp;


    if (left.empty()) {
        return vector<vector<int>>();
    }
    
    for (int i = 0; i < (m - k)/4; i++) {
        temp0 = getStretchedSolutionH(n, 4);
        if (temp0.empty()) {
            return vector<vector<int>>();
        }
        
        temp = join(temp, temp0, n, (4 * i), n, 4, horizontal);
        
        if (temp.empty()) {
            return vector<vector<int>>();
        }
    }
    
    return join(left, temp, n, k, n, (m - k), horizontal);
}

vector<vector<int>> thirdCaseSolution(int n, int m) {
    if (n != 4 || m <= 10) {
        return vector<vector<int>>();
    }
    
    int k = ((m - 6) % 5) + 6;
    
    vector<vector<int>> left = getStretchedSolutionV(n, k);
    
    vector<vector<int>> temp0;
    
    vector<vector<int>> temp;

    if (left.empty()) {
        return vector<vector<int>>();
    }
    
    for (int i = 0; i < (m - k)/5; i++) {
        temp0 = getDoubleLoopSolution(n, 5);
        
        if (temp0.empty()) {
            return vector<vector<int>>();
        }
        
        temp = join(temp, temp0, n, (5 * i), n, 5, horizontalDoubleLoop);
        
        //printBoard(temp);
        
        if (temp.empty()) {
            return vector<vector<int>>();
        }
    }
    
    return join(left, temp, n, k, n, (m - k), horizontalStretchedDoubleLoop);
}

vector<vector<int>> fourthCaseSolution1(int n, int m, stretchedType stretchedType) {
    int m1 = (m/4) * 2 + m % 2;
    int m2 = m - m1;
    
    vector<vector<int>> left = optimalOpenKnight(n, m1, stretchedType);
    vector<vector<int>> right = optimalOpenKnight(n, m2, horizontalStretched);
    
    if (left.empty() || right.empty()) {
        return vector<vector<int>>();
    }
    
    return join(left, right, n, m1, n, m2, horizontal);
}

vector<vector<int>> fourthCaseSolution2(int n, int m, stretchedType stretchedType) {
    int n1 = (n/4) * 2 + n % 2;
    int n2 = n - n1;
    
    vector<vector<int>> top = optimalOpenKnight(n1, m, stretchedType);
    vector<vector<int>> bottom = optimalOpenKnight(n2, m, verticalStretched);
    
    if (top.empty() || bottom.empty()) {
        return vector<vector<int>>();
    }
    
    return join(top, bottom, n1, m, n2, m, vertical);
}

vector<vector<int>> fifthCaseSolution(int n, int m, stretchedType stretchedType) {
    int n1 = (n/4) * 2 + n % 2;
    int n2 = n - n1;
    
    int m1 = (m/4) * 2 + m % 2;
    int m2 = m - m1;
    
    vector<vector<int>> topLeft = optimalOpenKnight(n1, m1, anyStretched);
    vector<vector<int>> topRight = optimalOpenKnight(n1, m2, horizontalStretched);
    vector<vector<int>> bottomLeft = optimalOpenKnight(n2, m1, verticalStretched);
    vector<vector<int>> bottomRight = optimalOpenKnight(n2, m2, horizontalStretched);
    
    if (topLeft.empty() || topRight.empty() || bottomLeft.empty() || bottomRight.empty()) {
        return vector<vector<int>>();
    }
    
    vector<vector<int>> top = join(topLeft, topRight, n1, m1, n1, m2, horizontal);
    vector<vector<int>> bottom = join(bottomLeft, bottomRight, n2, m1, n2, m2, horizontal);
    
    return join(top, bottom, n1, m, n2, m, vertical);
}

vector<vector<int>> optimalOpenKnight(int N, int M, stretchedType stretchedType) {
    int n = N;
    int m = M;
        
    if (n <= 10 && m <= 10) {
        return firstCaseSolution(n, m, stretchedType);
    }
    else if (n == 3 && m > 10) {
        return secondCaseSolution(n, m);
    }
    else if (n == 4 && m > 10) {
        return thirdCaseSolution(n, m);
    }
    else if (n >= 5 && n <= 10 && m > 10) {
        return fourthCaseSolution1(n, m, stretchedType);
    }
    //else if (m >= 5 && m <= 10 && n > 10) {
    //    return fourthCaseSolution2(n, m, stretchedType);
    //}
    else {
        return fifthCaseSolution(n, m, stretchedType);
    }
}

int main() {
    initBaseSolutionsDoubleLoop();
    initBaseSolutionsOpen();
    initBaseSolutionsStretchedV();
    initBaseSolutionsStretchedH();
    
    int N, M;
    bool needPrintBoard;

    cout << "Print board? (1 or 0):" << endl;
    cin >> needPrintBoard;
    
    cout << "Enter the dimensions of the board (N M): " << endl;
    cin >> N >> M;
    
    auto start = chrono::high_resolution_clock::now();
    
    vector<vector<int>> board = optimalOpenKnight(min(N,M), max(N,M), noStretched);
    
    auto end = chrono::high_resolution_clock::now();
    
    chrono::duration<double> duration = end - start;
    
    if (board.empty()) {
        cout << "No solution exists for the given board size.\n";
    }
    else {
        cout << "Knight's tour completed successfully:\n";
        if (needPrintBoard) {
            printBoard(board);
        }
        else {
            cout << "Print board is disabled" << endl;
        }
    }
    
    cout << "\n\nExecution time: " << duration.count() << " seconds" << endl;

    return 0;
}

