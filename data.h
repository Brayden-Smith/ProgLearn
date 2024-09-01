#include "LA.h"
#include <opencv2/opencv.hpp>
#include <fstream>

Matrix flattenMatrix(Matrix& matrixToFlatten);

Matrix shuffleRowOrder(Matrix& matrix);

std::vector<Matrix> trainTestSplit(Matrix& data, double percent); // Data should be in rows, percent is the percent (0.0 - 1.0) of data that should be in the training set

Matrix convertCVMatrixToMatrix(const cv::Mat& CVMatrix); // Converts a single-channel OpenCV matrix to a matrix object

Matrix importGrayscaleImage(const std::string& path); // Returns a grayscale image from a given path as an internal matrix object

std::vector<Matrix> importGrayscaleImageFamily(const std::string& path); /* Returns a vector of grayscale images found in the same file. The "path" argument is the path to a text file
containing a list of image file paths */

Matrix vectorToMatrixOfMatrices(std::vector<Matrix>& vectorOfMatrices); // Returns a matrix (order-preserving) with rows consisting of flattened matrices from the given vector





