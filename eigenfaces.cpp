#include "eigenfaces.h"
#include <iostream>

Matrix meanFlattenedFace(Matrix matrixOfFlattenedFaces) {
    Matrix matrixToReturn(1, matrixOfFlattenedFaces.getNumCols());

    for (int i = 0; i < matrixOfFlattenedFaces.getNumCols(); i++) {
        ComplexNum sum(0, 0);
        for (int j = 0; j < matrixOfFlattenedFaces.getNumRows(); j++) {
            sum += matrixOfFlattenedFaces(j, i);
        }
        matrixToReturn(0, i) = sum / matrixOfFlattenedFaces.getNumRows();
    }

    return matrixToReturn;
}


FaceSpace::FaceSpace(Matrix faceMatrixToAnalyze, double percentVariance) : faceDatabase(faceMatrixToAnalyze), eigenFaceDatabase(1, 1), numFaces(faceMatrixToAnalyze.getNumRows()), numEigenFaces(0), numFacesAdded(0), numFacesRemoved(0), numOfTopEigenFaces(0) {
    ComplexNum percentVarianceCNum(percentVariance, 0);
    // Must now create the eigenfaces
    // Normalize data
    Matrix faceMatrix = faceMatrixToAnalyze;
    ComplexNum number(1.0/255, 0);

    faceMatrix = faceMatrix * number;
    //std::cout << "Normed face matrix:\n" << faceMatrix << std::endl;

    meanFace = meanFlattenedFace(faceMatrix);

    //Matrix newMatrix = meanFace;
    //std::cout << "Mean flattened face:\n" << newMatrix << std::endl;

    // Subtract mean face from faceMatrix
    Matrix meanSubtractedMatrix = faceMatrix;
    for (int i = 0; i < faceMatrix.getNumRows(); i++) {
        for (int j = 0; j < faceMatrix.getNumCols(); j++) {
            meanSubtractedMatrix(i, j) = meanSubtractedMatrix(i, j) - meanFace(0, j);
        }
    }



    Matrix littleCovarianceMatrix = littleCovariance(meanSubtractedMatrix);

    //std::cout << "Little covariance matrix:\n" << littleCovarianceMatrix << std::endl;


    // Compute small eigenfaces
    std::vector<ComplexNum> eigenvalues;

    std::cout << "We try to compute eigenvectors" << std::endl;
    std::cout << "little covariance matrix has dimensions " << littleCovarianceMatrix.getNumRows() << " x " << littleCovarianceMatrix.getNumCols() << std::endl;

    std::vector<Matrix> smallEigenFaces = eigenvectors(&littleCovarianceMatrix, eigenvalues);
    std::cout << "We computed eigenvectors" << std::endl;

    for (int i = 0; i < eigenvalues.size(); i++) {
        std::cout << eigenvalues[i] << std::endl;
    }

    // Obtain final eigenfaces via linear transformation
    std::vector<Matrix> eigenFaces; // Column vectors
    eigenFaces.reserve(smallEigenFaces.size());

    Matrix meanMatrixTranspose = conjTranspose(&meanSubtractedMatrix);
    for (int i = 0; i < smallEigenFaces.size(); i++) {
        Matrix eigenFace = matMul(&meanMatrixTranspose, &smallEigenFaces[i]); // Double-check since there might be some reference weirdness in the second argument

        // Optional, but strongly recommended
        normalizeVectorsInMatrix(&eigenFace);
        eigenFaces.push_back(eigenFace);
    }

    // Select the top k eigenfaces using the given percentVariance argument

    std::reverse(eigenvalues.begin(), eigenvalues.end());
    std::reverse(eigenFaces.begin(), eigenFaces.end());

    ComplexNum totalVariance(0, 0);
    for (int i = 0; i < eigenvalues.size(); i++) {
        totalVariance += eigenvalues[i];
    }
    std::cout << "Total variance is: " << totalVariance << std::endl;

    ComplexNum eigenSum(0, 0);
    int k = 0;

    while (eigenSum/totalVariance < percentVarianceCNum) {
        eigenSum += eigenvalues[k];
        k += 1;
    }

    // Build face matrix of eigenfaces
    Matrix eigenFaceDatabaseBuilder(eigenFaces.size(), meanMatrixTranspose.getNumRows()); // Double-check dimensions are correct

    for (int i = 0; i < eigenFaces.size(); i++) {
        Matrix rowFace = flattenMatrix(eigenFaces[i]);
        eigenFaceDatabaseBuilder(i) = rowFace;
    }

    // Fill out the object information
    eigenFaceDatabase = eigenFaceDatabaseBuilder;
    numEigenFaces = eigenFaceDatabase.getNumRows();
    numOfTopEigenFaces = k;
}



FaceSpace::FaceSpace(std::string& path, double percentVariance) : faceDatabase(Matrix(1, 1)), meanFace(Matrix(1,1)), eigenFaceDatabase(Matrix(1, 1)), numFaces(0), numEigenFaces(0), numFacesAdded(0), numFacesRemoved(0), numOfTopEigenFaces(0) {
    // todo
}




long FaceSpace::getNumFaces() {
    return numFaces;
}

unsigned int FaceSpace::getNumEigenFaces() {
    return numEigenFaces;
}

unsigned int FaceSpace::getNumTopEigenFaces() {
    return numOfTopEigenFaces;
}

unsigned int FaceSpace::getNumFacesAdded() {
    return numFacesAdded;
}

unsigned int FaceSpace::getNumFacesRemoved() {
    return numFacesRemoved;
}

Matrix FaceSpace::getFaceEntry(int n) {
    if (numFaces < n) {
        throw std::invalid_argument("FaceSpace::getFaceEntry: argument out of bounds!");
    }
    return faceDatabase[n - 1];
}

Matrix FaceSpace::getMeanFace() {
    return meanFace;
}

Matrix FaceSpace::getEigenFaceEntry(int n) {
    if (numEigenFaces < n) {
        throw std::invalid_argument("FaceSpace::getEigenFaceEntry: argument out of bounds!");
    }
    return eigenFaceDatabase[n - 1];
}







