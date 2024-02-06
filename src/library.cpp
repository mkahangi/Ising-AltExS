#include <vector>
#include <stdexcept>

#include "function.hpp"

//converting radians to degrees
double radiansToDegrees(double radians) {
    return radians * (180.0 / M_PI);
}

//******************************************************************************************************************//

// Function to convert function of T  to a vector and store it.
void storevector(vector<double>&vec,int row, double rowstep, double (*function)(double)) {
    for (int i = 0; i < row; i+=1) {
        double Tval = static_cast<double>(i) * rowstep;
        vec[i] = function(Tval);
        }
}

// Function to convert a function of T and muB to a 2D matrix and store it.
void storeMatrix(vector<vector<double>>& matrix, int row, double rowstep, int col, double colstep, double (*function)(double, double)) {
    clock_t tStart = clock(); // Timer start
    // Create a result matrix with the same dimensions as the input matrix
    for (int i = 0; i < row; i += 1) {
        for (int j = 0; j < col; j += 1) {
            double Tval = static_cast<double>(i) * rowstep;
            double muBval = static_cast<double>(j) * colstep;

            matrix[i][j] = function(Tval, muBval);
           
            // Calculate progress percentage
            double progress = static_cast<double>(i * col + j) / static_cast<double>(row * col) * 100;

            // Print the progress bar
            cout << "\rProgress: [" << setw(3) << static_cast<int>(progress) << "%] [";
            int barWidth = 50;
            int progressWidth = static_cast<int>(progress / 100 * barWidth);
            for (int k = 0; k < barWidth; ++k) {
                cout << (k < progressWidth ? '=' : ' ');
            }
            cout << "]" << flush;
        }
    }
    printf(" Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    // Print a newline to complete the progress bar
    cout << endl;
}



std::vector<std::vector<double>> integrateMatrix(const std::vector<std::vector<double>>& matrix, int row, double rowstep, int col, double colstep, bool direction) {
    // Parameter validation
    if (matrix.empty() || matrix.size() != row || matrix[0].size() != col) {
        throw std::invalid_argument("Invalid matrix dimensions");
    }

    std::vector<std::vector<double>> result(row, std::vector<double>(col, 0.0));

    // Helper lambda for trapezoidal rule integration
    auto trapezoidalRule = [](double prevResult, double step, double current, double previous) {
        // Check if either current or previous is NaN and handle accordingly
        if (std::isnan(current) || std::isnan(previous)) {
            if (std::isnan(current) && std::isnan(previous)) {
                // Both are NaN, skip this step entirely
                return prevResult;
            } else if (std::isnan(current)) {
                // Only current is NaN, use only the previous value
                return prevResult + step * previous; // This might need to be adjusted based on how you want to handle a single value
            } else {
                // Only previous is NaN, start fresh from current
                return prevResult + step * current; // Similarly, might need adjustment
            }
        }
        return prevResult + 0.5 * step * (current + previous);
    };

    if (direction) {  // Row-wise integration
        for (int i = 0; i < row; ++i) {
            if (!std::isnan(matrix[i][0])) {
                result[i][0] = matrix[i][0] * rowstep;
            }
            for (int j = 1; j < col; ++j) {
                if (!std::isnan(matrix[i][j]) || !std::isnan(matrix[i][j - 1])) {
                    result[i][j] = trapezoidalRule(result[i][j - 1], colstep, matrix[i][j], matrix[i][j - 1]);
                }
                // If both current and previous are NaN, result[i][j] remains as initialized to 0.0
            }
        }
    } else {  // Column-wise integration
        for (int j = 0; j < col; ++j) {
            if (!std::isnan(matrix[0][j])) {
                result[0][j] = matrix[0][j] * colstep;
            }
            for (int i = 1; i < row; ++i) {
                if (!std::isnan(matrix[i][j]) || !std::isnan(matrix[i - 1][j])) {
                    result[i][j] = trapezoidalRule(result[i - 1][j], rowstep, matrix[i][j], matrix[i - 1][j]);
                }
                // If both are NaN, then result[i][j] remains as initialized
            }
        }
    }

    return result;
}


// std::vector<std::vector<double>> integrateMatrix(const std::vector<std::vector<double>>& matrix, int row, double rowstep, int col, double colstep, bool direction) {
//     // Parameter validation
//     if (matrix.empty() || matrix.size() != row || matrix[0].size() != col) {
//         throw std::invalid_argument("Invalid matrix dimensions");
//     }

//     std::vector<std::vector<double>> result(row, std::vector<double>(col, 0.0));

//     // Helper lambda for trapezoidal rule integration
//     auto trapezoidalRule = [](double prevResult, double step, double current, double previous) {
//         return prevResult + 0.5 * step * (current + previous);
//     };

//     if (direction) {  // Row-wise integration
//         for (int i = 0; i < row; ++i) {
//             result[i][0] = matrix[i][0] * rowstep;
//             for (int j = 1; j < col; ++j) {
//                 result[i][j] = trapezoidalRule(result[i][j - 1], colstep, matrix[i][j], matrix[i][j - 1]);
//             }
//         }
//     } else {  // Column-wise integration
//         for (int j = 0; j < col; ++j) {
//             result[0][j] = matrix[0][j] * colstep;
//             for (int i = 1; i < row; ++i) {
//                 result[i][j] = trapezoidalRule(result[i - 1][j], rowstep, matrix[i][j], matrix[i - 1][j]);
//             }
//         }
//     }

//     return result;
// }



// // Function to compute the derivative of a 2D matrix along rows or columns
vector<vector<double>> deriv_matrix(const vector<vector<double>>& matrix, int rowSize, double rowStep, int colSize, double colStep, int derivativeDirection) {
    // Create a result matrix with the same dimensions as the input matrix
    vector<vector<double>> result(rowSize, vector<double>(colSize, 0.0));

    if (derivativeDirection == 1) {  // Derivative along rows
        for (int i = 0; i < rowSize; ++i) {
            for (int j = 0; j < colSize; ++j) {
                // Compute the derivative along rows using central difference for interior points
                result[i][j] = (j > 0 && j < colSize - 1) ? (matrix[i][j + 1] - matrix[i][j - 1]) / (2 * colStep) :
                                                           (j == 0) ? (matrix[i][1] - matrix[i][0]) / colStep :
                                                                      (matrix[i][j] - matrix[i][j - 1]) / colStep;
            }
        }
    } else {  // Derivative along columns
        for (int i = 0; i < rowSize; ++i) {
            for (int j = 0; j < colSize; ++j) {
                // Compute the derivative along columns using central difference for interior points
                result[i][j] = (i > 0 && i < rowSize - 1) ? (matrix[i + 1][j] - matrix[i - 1][j]) / (2 * rowStep) :
                                                           (i == 0) ? (matrix[1][j] - matrix[0][j]) / rowStep :
                                                                      (matrix[i][j] - matrix[i - 1][j]) / rowStep;
            }
        }
    }

    return result;
}




