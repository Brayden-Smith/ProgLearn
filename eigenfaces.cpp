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

    Matrix faceMatrix = faceMatrixToAnalyze;
    ComplexNum number(1.0/255, 0);

    faceMatrix = faceMatrix * number;

    meanFace = meanFlattenedFace(faceMatrix);

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

    std::vector<Matrix> smallEigenFaces = CVEigenvectors(&littleCovarianceMatrix, eigenvalues);
    //std::vector<Matrix> smallEigenFaces = eigenvectors(&littleCovarianceMatrix, eigenvalues);
    std::cout << "We computed eigenvectors" << std::endl;


    for (int i = 0; i < eigenvalues.size(); i++) {
        std::cout << eigenvalues[i] << std::endl;
    }

    std::cout << "Num eigenvectors is " << smallEigenFaces.size() << std::endl;

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

    for (int i = 0; i < eigenvalues.size(); i++) {
        std::cout << eigenvalues[i] << std::endl;
    }

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

    std::cout << "Face builder dims: " << eigenFaces.size() << "x" << meanMatrixTranspose.getNumRows() << std::endl;

    for (int i = 0; i < eigenFaces.size(); i++) {
        Matrix rowFace = flattenMatrix(eigenFaces[i]);
        for (int j = 0; j < rowFace.getNumCols(); j++) {
            eigenFaceDatabaseBuilder(i, j) = rowFace(0, j);
        }
        //eigenFaceDatabaseBuilder(i) = rowFace;
    }



    // Fill out the object information
    eigenFaceDatabase = eigenFaceDatabaseBuilder;
    numEigenFaces = eigenFaceDatabase.getNumRows();
    numOfTopEigenFaces = k;
    std::cout << "Done with face object construction" << std::endl;
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

    Matrix matrixToReturn(1, faceDatabase.getNumCols());

    for (int i = 0; i < faceDatabase.getNumCols(); i++) {
        matrixToReturn(0, i) = faceDatabase(n-1, i);
    }

    return matrixToReturn;
}

Matrix FaceSpace::getMeanFace() {
    return meanFace;
}

Matrix FaceSpace::getEigenFaceEntry(int n) {
    if (numEigenFaces < n) {
        throw std::invalid_argument("FaceSpace::getEigenFaceEntry: argument out of bounds!");
    }
    Matrix matrixToReturn(1, eigenFaceDatabase.getNumCols());

    for (int i = 0; i < eigenFaceDatabase.getNumCols(); i++) {
        matrixToReturn(0, i) = eigenFaceDatabase(n-1, i);
    }

    return matrixToReturn;
}


Matrix FaceSpace::matchFace(Matrix &faceToMatch, double similarityScore) {
    Matrix faceToMatchColumnVector = conjTranspose(&faceToMatch);
    Matrix columnMeanFace = conjTranspose(&meanFace);
    Matrix processedFaceColumnVector = faceToMatchColumnVector - columnMeanFace;

    Matrix faceToMatchWeightVector(numEigenFaces, 1);
    for (int i = 0; i < numEigenFaces; i++) {
        Matrix flattenedEigenFace = eigenFaceDatabase[i];
        Matrix columnEigenFace = conjTranspose(&flattenedEigenFace);
        ComplexNum projectionResult = innerProduct(columnEigenFace, processedFaceColumnVector);
        faceToMatchWeightVector(i, 0) = projectionResult;
    }

    normalizeVectorsInMatrix(&faceToMatchWeightVector);

    // Now go through each vector, and compare

    double minDistance = 2.1;
    int minDistFaceIndex = 0;
    for (int i = 0; i < numFaces; i++) {
        // todo


    }



    return Matrix();
}







