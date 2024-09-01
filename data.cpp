#include "data.h"
#include <random>
#include <opencv2/opencv.hpp>
#include <fstream>


Matrix flattenMatrix(Matrix& matrixToFlatten) {
    Matrix flattenedMatrix(matrixToFlatten.getNumRows() * matrixToFlatten.getNumCols(), 1);
    int cumulativeSum = 0;
    for (int i = 0; i < matrixToFlatten.getNumRows(); i++) {
        for (int j = 0; j < matrixToFlatten.getNumCols(); j++) {
            flattenedMatrix(cumulativeSum++, 0) = matrixToFlatten(i, j);
        }
    }
    return flattenedMatrix;
}


Matrix shuffleRowOrder(Matrix& matrix) {
    return Matrix(3, 3);
} // todo


std::vector<Matrix> trainTestSplit(Matrix& data, double percent) {
    if (percent > 1 || percent < 0) {
        throw std::invalid_argument("trainTestSplit: percent must be between 0 and 1!");
    }

    int numRows = data.getNumRows();
    int numTrainRows = std::floor(percent * numRows);
    int numTestRows = numRows - numTrainRows;

    Matrix trainMatrix(numTrainRows, data.getNumCols());
    Matrix testMatrix(numTestRows, data.getNumCols());


    for (int i = 0; i < numTrainRows; i++) {
        for (int j = 0; j < data.getNumCols(); j++) {
            trainMatrix(i, j) = data(i, j);
        }
    }

    for (int i = numTrainRows; i < numRows; i++) {
        for (int j = 0; j < data.getNumCols(); j++) {
            testMatrix(i - numTrainRows, j) = data(i, j);
        }
    }

    return {trainMatrix, testMatrix};
}



Matrix convertCVMatrixToMatrix(const cv::Mat& CVMatrix) { // Converts a single-channel OpenCV matrix to a matrix object
    if (CVMatrix.channels() != 1) {
        throw std::invalid_argument("convertCVMatrixToMatrix: CVMatrix argument must only have one channel!");
    }
    if (CVMatrix.empty()) {
        throw std::invalid_argument("convertCVMatrixToMatrix: Matrix is empty!");
    }
    if (CVMatrix.type() != CV_8U) {
        throw std::invalid_argument("CVMatrix must be of type CV_8U (0-255)!");
    }

    int numRows = CVMatrix.rows;
    int numCols = CVMatrix.cols;
    Matrix matrixToReturn(numRows, numCols);

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            double valueInCVMatrix = CVMatrix.at<uchar>(i, j);
            ComplexNum CNumValue(valueInCVMatrix, 0);
            matrixToReturn(i, j) = CNumValue;
        }
    }

    return matrixToReturn;
}


Matrix importGrayscaleImage(const std::string& path) {
    cv::Mat CVImage = cv::imread(path, cv::IMREAD_GRAYSCALE); // Read the image from the given path in grayscale and return a cv matrix

    if (CVImage.empty()) {
        throw std::runtime_error("Failed to load image from: " + path);
    }
    if (CVImage.channels() != 1) {
        throw std::runtime_error("Image is not grayscale: " + path);
    }
    if (CVImage.type() != CV_8U) { // If grayscale image is not in 8-bit unsigned integer format, convert it to this format
        double min, max;
        cv::minMaxLoc(CVImage, &min, &max);
        CVImage.convertTo(CVImage, CV_8U, 255.0 / (max - min), (-1 * min * 255.0)/(max - min));
    }

    return convertCVMatrixToMatrix(CVImage);

}


std::vector<Matrix> importGrayscaleImageFamily(const std::string& path) {
    std::vector<Matrix> vectorToReturn;
    std::ifstream file(path);
    std::string imagePath;

    if (!file.is_open()) {
        throw std::runtime_error("importGrayscaleImage: Could not open file: " + path);
    }

    while (std::getline(file, imagePath)) {
        try {
            Matrix imageMatrixToAppend = importGrayscaleImage(imagePath);
            vectorToReturn.push_back(imageMatrixToAppend);
        } catch (const std::exception& error) {
            std::cout << "importGrayscaleImage: Error processing image " << imagePath << ": " << error.what() << std::endl;
        }

    }
    file.close();

    if (vectorToReturn.empty()) {
        throw std::runtime_error("importGrayscaleImage: No valid images could be imported from " + imagePath);
    }

    return vectorToReturn;
}



Matrix vectorToMatrixOfMatrices(std::vector<Matrix>& vectorOfMatrices) {
    if (vectorOfMatrices.empty()) {
        throw std::invalid_argument("vectorToMatrixOfMatrices: Vector does not contain any matrices!");
    }

    Matrix matrixToReturn(vectorOfMatrices.size(), vectorOfMatrices[0].getNumRows() * vectorOfMatrices[0].getNumCols());
    int numCols = matrixToReturn.getNumCols();

    for (int i = 0; i < vectorOfMatrices.size(); i++) { // Will need to change to long type if dealing with very many large image files
        if (vectorOfMatrices[i].getNumRows() * vectorOfMatrices[i].getNumCols() != numCols) {
            throw std::invalid_argument("vectorToMatrixOfMatrices: Matrices in vector not of the same size!");
        }

        Matrix flattenedMatrix = flattenMatrix(vectorOfMatrices[i]);
        matrixToReturn(i) = flattenedMatrix; // Set row i of matrixToReturn equal to the flattened matrix in position i of the vector of matrices
    }

    return matrixToReturn;
}


