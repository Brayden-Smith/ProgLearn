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



FaceSpace::FaceSpace(Matrix faceMatrix, double percentVariance) : faceDatabase(faceMatrix), meanFace(meanFlattenedFace(faceMatrix)), eigenFaceDatabase(1, 1), numFaces(faceMatrix.getNumRows()), numEigenFaces(0), numFacesAdded(0), numFacesRemoved(0), numOfTopEigenFaces(0) {
    ComplexNum percentVarianceCNum(percentVariance, 0);
    // Must now create the eigenfaces

    // Subtract mean face from faceMatrix
    Matrix meanSubtractedMatrix = faceMatrix;
    for (int i = 0; i < faceMatrix.getNumRows(); i++) {
        for (int j = 0; j < faceMatrix.getNumCols(); j++) {
            meanSubtractedMatrix(i, j) = meanSubtractedMatrix(i, j) - meanFace(0, j);
        }
    }

    Matrix littleCovarianceMatrix = littleCovariance(meanSubtractedMatrix);

    // Compute small eigenfaces
    std::vector<ComplexNum> eigenvalues;
    std::vector<Matrix> smallEigenFaces = eigenvectors(&littleCovarianceMatrix, eigenvalues);

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

    ComplexNum totalVariance(0, 0);
    for (int i = 0; i < eigenvalues.size(); i++) {
        totalVariance += eigenvalues[i];
    }

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




