#include<algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <numeric>

// import itertools
// import numpy as np

// # Edge ranking method
// # We give each edge a ranking and normalise from best to worst scenario.

// def calculate_LH(x_L, x_H, weights):
//     """Given the best and worst case scenarios, x_L and x_H, we calculate the unit weight for the CPT."""
//     return x_L/sum(weights), x_H/sum(weights)

float sumOfVector(std::vector<float> weights){
    return std::accumulate(weights.begin(), weights.end(), 0.0f);
}

float divideBySumOfVector(float number, std::vector<float> weights) {
    // Create a vector of floats
    // std::vector<float> myVector = {1.5, 2.3, 4.7, 5.1, 0.9};

    return number / sumOfVector(weights);
}

// def combinations(num_weights):
//     return list(itertools.product([False, True], repeat=num_weights))

std::vector<std::vector<int>> generateBinaryCombinations(int length) {
    // Calculate the total number of combinations (2^length)
    int numCombinations = pow(2, length);

    // Create the matrix to store the combinations
    std::vector<std::vector<int>> combinations(numCombinations, std::vector<int>(length));

    // Iterate through each combination
    for (int i = 0; i < numCombinations; ++i) {
        // Generate the binary representation of the current combination index
        for (int j = 0; j < length; ++j) {
            combinations[i][j] = (i >> j) & 1;
        }
    }

    return combinations;
}

// def create_boolean_matrix(x_L, x_H, num_weights):
//     # Generate all possible combinations of True and False for the given number
//     all_combinations = combinations(num_weights)
//     # Create a boolean matrix with each combination as a row
//     boolean_matrix = [[x_L if bit else x_H for bit in combination] for combination in all_combinations]

//     return np.array(boolean_matrix)


std::vector<std::vector<float>> replaceValuesInMatrix(
    std::vector<std::vector<int>> intMatrix
    , float low
    , float high) {
    
    std::vector<std::vector<float>> floatMatrix;
    
    for (const std::vector<int> intRow : intMatrix) {
        // Create a new row for the floatMatrix
        std::vector<float> floatRow;

        // Iterate through each element in the intRow and convert to float
        for (int value : intRow) {
            if (value == 0) {
                floatRow.push_back(high);
            } else if (value == 1) {
                floatRow.push_back(low);
            }
            ;
        }

        // Add the floatRow to the floatMatrix
        floatMatrix.push_back(floatRow);
    }
    return floatMatrix;
}

// def calculate_CPT(weights, x_L, x_H):
//     boolean_matrix = create_boolean_matrix(
//         *calculate_LH(x_L, x_H, weights),
//         len(weights)
//     )

//     return sum(np.transpose(np.array([weights]*len(boolean_matrix))*boolean_matrix))
std::vector<float> calculateCPT(std::vector<float> weights, float low, float high){
    
    std::vector<float> cpt;

    std::vector<std::vector<float>> matrix = replaceValuesInMatrix(
        generateBinaryCombinations(weights.size())
        , divideBySumOfVector(low, weights)
        , divideBySumOfVector(high, weights)
    );

    std::size_t rows = matrix.size();
    std::size_t cols = matrix[0].size();  // Assuming all rows have the same size

    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            matrix[i][j] = matrix[i][j] * weights[j];  // Multiply each element with itself
        }
        cpt.push_back(sumOfVector(matrix[i]));
    }
    return cpt;
}

// if __name__ == "__main__":
//     weights = [1, 2]
//     x_L = 0.1
//     x_H = 0.9
//     matrix = calculate_CPT(weights, x_L, x_H)
//     np.savetxt("low", matrix, newline='\t', fmt="%1f")
//     np.savetxt("high", 1-matrix, newline='\t', fmt="%1f")

int main(){
    // Print the result
    std::vector<float> testVector = {2, 3, 5};
    std::cout << "test vector is {";
    for (int i: testVector)
    {
        std::cout << i << ", ";
    };
    std::cout << "}\n";
    std::cout << "The total sum of the vector is: " << 1/divideBySumOfVector(1, testVector) << std::endl;
    std::cout << "\n";
    std::cout << "if we set the number to 2...\n" << divideBySumOfVector(2, testVector) << std::endl;
    std::cout << "\n";

    std::cout << "For 3 entires there should be 8 combinations (2^3):\n";

    std::vector<std::vector<int>> combination = generateBinaryCombinations(3);
    int line = 0;

    for (std::vector<int> i: combination)
    {
        line += 1;
        std::cout << line <<": [ ";
        for (int j: i) 
        {
            std::cout << "" << j << ", ";
        }
        std::cout << "]\n";
    }

    std::vector<std::vector<float>> matrix = replaceValuesInMatrix(combination, divideBySumOfVector(0.9, testVector), divideBySumOfVector(0.1, testVector));

    std::size_t rows = matrix.size();
    std::size_t cols = matrix[0].size();  // Assuming all rows have the same size

    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            matrix[i][j] = matrix[i][j] * testVector[j];  // Multiply each element with itself
        }
    }

    float rowSum;
    line = 0;
    std::cout << "\n";
    for (std::vector<float> i: matrix)
    {
        line += 1;
        std::cout << line << ": [ ";
        for (float j: i) 
        {
            std::cout << "" << j << ", ";
        }
        rowSum = sumOfVector(i);
        std::cout << "] " << rowSum << " | " << 1-rowSum << "\n";
    }

    for (float i: calculateCPT(testVector, 0.1, 0.9))
    {
        std::cout << i << ", ";
    };
    std::cout << "\n";

    return 0;
}