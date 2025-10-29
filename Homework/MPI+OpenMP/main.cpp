//
//  main.cpp
//  KnightTour MPI
//
//  Created by Serguei Diaz on 21.12.2024.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <mpi.h>
#include <omp.h>


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
    horizontalStretched = 0,
    verticalStretched = 1,
    anyStretched = 2,
    noStretched = 3
};

vector<vector<int>> optimalOpenKnight(int N, int M, stretchedType stretchedType, int depth, int num_threads);
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
        
        if (temp.empty()) {
            return vector<vector<int>>();
        }
    }
    
    return join(left, temp, n, k, n, (m - k), horizontalStretchedDoubleLoop);
}

vector<vector<int>> fourthCaseSolution1(int n, int m, stretchedType stretchedType) {
    int m1 = (m/4) * 2 + m % 2;
    int m2 = m - m1;
    
    vector<vector<int>> left = optimalOpenKnight(n, m1, stretchedType, 0, 1);
    vector<vector<int>> right = optimalOpenKnight(n, m2, horizontalStretched, 0, 1);
    
    if (left.empty() || right.empty()) {
        return vector<vector<int>>();
    }
    
    return join(left, right, n, m1, n, m2, horizontal);
}

vector<vector<int>> fourthCaseSolution2(int n, int m, stretchedType stretchedType) {
    int n1 = (n/4) * 2 + n % 2;
    int n2 = n - n1;
    
    vector<vector<int>> top = optimalOpenKnight(n1, m, stretchedType, 0, 1);
    vector<vector<int>> bottom = optimalOpenKnight(n2, m, verticalStretched, 0, 1);
    
    if (top.empty() || bottom.empty()) {
        return vector<vector<int>>();
    }
    
    return join(top, bottom, n1, m, n2, m, vertical);
}

vector<vector<int>> fifthCaseSolution(int n, int m, stretchedType stretchedType, int depth, int num_threads) {
    int n1 = (n/4) * 2 + n % 2;
    int n2 = n - n1;
    
    int m1 = (m/4) * 2 + m % 2;
    int m2 = m - m1;
    
    vector<vector<int>> topLeft, topRight, bottomLeft, bottomRight;
    
    if (((depth - 1) * 4) + 1 < num_threads && depth != 0) {
        #pragma omp parallel
        {
            #pragma omp single
            {
                #pragma omp task
                topLeft = optimalOpenKnight(n1, m1, anyStretched, depth + 1, num_threads);

                #pragma omp task
                topRight = optimalOpenKnight(n1, m2, horizontalStretched, depth + 1, num_threads);

                #pragma omp task
                bottomLeft = optimalOpenKnight(n2, m1, verticalStretched, depth + 1, num_threads);

                #pragma omp task
                bottomRight = optimalOpenKnight(n2, m2, horizontalStretched, depth + 1, num_threads);

                #pragma omp taskwait
            }
        }
        
    }
    else {
        topLeft = optimalOpenKnight(n1, m1, anyStretched, depth, num_threads);
        topRight = optimalOpenKnight(n1, m2, horizontalStretched, depth, num_threads);
        bottomLeft = optimalOpenKnight(n2, m1, verticalStretched, depth, num_threads);
        bottomRight = optimalOpenKnight(n2, m2, horizontalStretched, depth, num_threads);
    }
    
    if (topLeft.empty() || topRight.empty() || bottomLeft.empty() || bottomRight.empty()) {
        return vector<vector<int>>();
    }
    
    vector<vector<int>> top = join(topLeft, topRight, n1, m1, n1, m2, horizontal);
    vector<vector<int>> bottom = join(bottomLeft, bottomRight, n2, m1, n2, m2, horizontal);
    
    return join(top, bottom, n1, m, n2, m, vertical);
}

vector<vector<int>> optimalOpenKnight(int N, int M, stretchedType stretchedType, int depth, int num_threads) {
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
        return fifthCaseSolution(n, m, stretchedType, depth, num_threads);
    }
}

vector<vector<vector<int>>> preCalculationForMPI(int n, int m, int nodes) {
    vector<vector<int>> topLeft, topRight, bottomLeft, bottomRight;
    vector<vector<vector<int>>> response;
    vector<vector<int>> tempCuadrants(4, vector<int>(2));
    
    int currentParts = 4;
    int n1, n2, m1, m2, prevParts;
    
    if (n > 10 && m > 10) {
        n1 = (n/4) * 2 + n % 2;
        n2 = n - n1;
        m1 = (m/4) * 2 + m % 2;
        m2 = m - m1;
        
        topLeft.push_back({n1, m1});
        topRight.push_back({n1, m2});
        bottomLeft.push_back({n2, m1});
        bottomRight.push_back({n2, m2});
        response = {topLeft, topRight, bottomLeft, bottomRight};
    }
    else {
        return {};
    }
    
    while (currentParts < nodes) {
        prevParts = currentParts / 4;
        
        for (int i = 0; i < prevParts; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 2; ++k) {
                    if (response[j][i][k] <= 10) {
                        return response;
                    }
                }
            }
        }
        
        currentParts = currentParts * 4;
        
        for (int i = 0; i < prevParts; ++i) {
            for (int j = 0; j < 4; ++j) {
                n = response[j][0][0];
                m = response[j][0][1];
                
                n1 = (n/4) * 2 + n % 2;
                n2 = n - n1;
                m1 = (m/4) * 2 + m % 2;
                m2 = m - m1;
                
                response[0].push_back({n1, m1});
                response[1].push_back({n1, m2});
                response[2].push_back({n2, m1});
                response[3].push_back({n2, m2});
                
                response[j].erase(response[j].begin());
            }
        }
    }
    
    return response;
}

int main(int argc, char** argv) {
    int rank, size;
    int N;
    int M;
    int num_threads;
    bool needPrintBoard;
    
    chrono::time_point<chrono::high_resolution_clock> start, end;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    initBaseSolutionsDoubleLoop();
    initBaseSolutionsOpen();
    initBaseSolutionsStretchedV();
    initBaseSolutionsStretchedH();
    
    if (rank == 0) {
        cout << "Print board? (1 or 0):" << endl;
        cin >> needPrintBoard;
        cout << "Enter the number of threads: " << endl;
        cin >> num_threads;
        cout << "Enter the dimensions of the board (N M): " << endl;
        cin >> N >> M;
        start = chrono::high_resolution_clock::now();
    }
    
    int sendcounts[size];
    int displs[size];
    int vectorSize;
    vector<vector<vector<int>>> parts;
    
    if (rank == 0) {
        int n = min(N,M);
        int m = max(N,M);
        parts = preCalculationForMPI(n, m, size);
        vectorSize = parts[0].size();
        
        int totalParts = vectorSize * 4;
        int partsPerNode = totalParts / size;
        int residualParts = totalParts % size;
        int offset = 0;
        
        for (int i = 0; i < size; ++i) {
            sendcounts[i] = (partsPerNode + (i < residualParts ? 1 : 0)) * 4;
            displs[i] = offset;
            offset += sendcounts[i];
        }
    }
    
    MPI_Bcast(sendcounts, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_threads, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    cout << "rank " << rank <<" num_threads: " << num_threads << endl;
    
    omp_set_num_threads(num_threads);
    omp_set_max_active_levels(num_threads / 4 + 1);
    
    int rowsToReceive = sendcounts[rank] / 4;
    
    vector<int> flatMatrix;
    
    if (rank == 0) {
        for (int i = 0; i < vectorSize; ++i) {
            for (int j = 0; j < 4; ++j) {
                flatMatrix.push_back(parts[j][i][0]);
                flatMatrix.push_back(parts[j][i][1]);
                flatMatrix.push_back(j);
                flatMatrix.push_back(i);
            }
        }
    }
    
    vector<int> recvBuffer(rowsToReceive * 4);
    
    MPI_Scatterv(flatMatrix.data(), sendcounts, displs, MPI_INT, recvBuffer.data(), recvBuffer.size(), MPI_INT, 0, MPI_COMM_WORLD);
    
    vector<vector<int>> rowsReceived(rowsToReceive, vector<int>(4));
    for (int i = 0; i < rowsToReceive; ++i) {
        for (int j = 0; j < 4; ++j) {
            rowsReceived[i][j] = recvBuffer[i * 4 + j];
        }
    }
    
    vector<vector<vector<int>>> response;
    vector<vector<int>> temp;
    for (int i = 0; i < rowsToReceive; ++i) {
        
        if (rowsReceived[i][2] == 0 && rowsReceived[i][3] == 0) {
            temp = optimalOpenKnight(rowsReceived[i][0], rowsReceived[i][1], noStretched, 0, num_threads);
            if (temp.empty()) {
                cout << "rank " << rank << " : No solution" << endl;
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            else {
                response.push_back(temp);
            }
        }
        else {
            if (rowsReceived[i][2] == 0) {
                temp = optimalOpenKnight(rowsReceived[i][0], rowsReceived[i][1], anyStretched, 0, num_threads);
            }
            else if (rowsReceived[i][2] == 2) {
                temp = optimalOpenKnight(rowsReceived[i][0], rowsReceived[i][1], verticalStretched, 0, num_threads);
            }
            else {
                temp = optimalOpenKnight(rowsReceived[i][0], rowsReceived[i][1], horizontalStretched, 0, num_threads);
            }
            if (temp.empty()) {
                cout << "rank " << rank << " : No solution" << endl;
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            else {
                response.push_back(temp);
            }
        }
    }
    
    vector<int> flatResponse;
    
    for (int i = 0; i < rowsToReceive; ++i) {
        flatResponse.insert(flatResponse.end(), rowsReceived[i].begin(), rowsReceived[i].end());
        for (const auto& row : response[i]) {
            flatResponse.insert(flatResponse.end(), row.begin(), row.end());
        }
    }
    
    int recvcounts[size];
    int displsGather[size];
    int totalSize = 0;
    int tempRowsCounter = 0;
    vector<int> finalResult;
    
    if (rank == 0) {
        totalSize = 0;
        for (int i = 0; i < size; ++i) {
            recvcounts[i] = 0;
            for (int j = 0; j < sendcounts[i] / 4; ++j) {
                recvcounts[i] = recvcounts[i] + 4 + flatMatrix[tempRowsCounter*4] * flatMatrix[tempRowsCounter*4 + 1];
                tempRowsCounter = tempRowsCounter + 1;
            }
            displsGather[i] = totalSize;
            totalSize += recvcounts[i];
        }
        finalResult.resize(totalSize);
    }
    
    MPI_Gatherv(flatResponse.data(), flatResponse.size(), MPI_INT, finalResult.data(), recvcounts, displsGather, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        int fn, fm, fj, fi, skip, first;
        vector<vector<vector<vector<int>>>> resultParts = {{}, {}, {}, {}};
        skip = 0;
        for (int i = 0; i < parts[0].size(); ++i) {
            for (int j = 0; j < 4; ++j) {
                fn = finalResult[skip];
                fm = finalResult[skip + 1];
                fj = finalResult[skip + 2];
                fi = finalResult[skip + 3];
                
                resultParts[j].push_back({});
                
                for (int k = 0; k < fn; ++k) {
                    first = (skip + 4) + (k * fm);
                    vector<int> row(finalResult.begin() + first, finalResult.begin() + first + fm);
                    resultParts[j][i].push_back(row);
                }
                
                skip = skip + 4 + fn * fm;
            }
        }
        
        vector<vector<int>> topLeft, topRight, bottomLeft, bottomRight, top, bottom, complete;
        int counter = 0;
        int subPos;
        while (!resultParts[1].empty()) {
            topLeft = resultParts[0][0];
            topRight = resultParts[1][0];
            bottomLeft = resultParts[2][0];
            bottomRight = resultParts[3][0];
            
            resultParts[0].erase(resultParts[0].begin());
            resultParts[1].erase(resultParts[1].begin());
            resultParts[2].erase(resultParts[2].begin());
            resultParts[3].erase(resultParts[3].begin());
            
            
            subPos = counter % 4;
            
            top = join(topLeft, topRight, topLeft.size(), topLeft[0].size(), topRight.size(), topRight[0].size(), horizontal);
            bottom = join(bottomLeft, bottomRight, bottomLeft.size(), bottomLeft[0].size(), bottomRight.size(), bottomRight[0].size(), horizontal);
            
            complete = join(top, bottom, top.size(), top[0].size(), bottom.size(), bottom[0].size(), vertical);
            
            resultParts[subPos].push_back(complete);
            
            counter = counter + 1;
        }
        
        end = chrono::high_resolution_clock::now();
        
        chrono::duration<double> duration = end - start;
        
        if (resultParts[0][0].empty()) {
            cout << "No solution exists for the given board size.\n" << endl;
        }
        else {
            cout << "Knight's tour completed successfully:\n" << endl;
            if (needPrintBoard) {
                printBoard(resultParts[0][0]);
            }
            else {
                cout << "Print board is disabled" << endl;
            }
        }
        
        cout << "\n\nExecution time: " << duration.count() << " seconds" << endl;
        
    }
    
    MPI_Finalize();
    
    return 0;
}

