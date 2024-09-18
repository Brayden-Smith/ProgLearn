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




FaceSpace::FaceSpace(Matrix& faceMatrixToAnalyze, double percentVariance) :  eigenFaceDatabase(1, 1), numFaces(faceMatrixToAnalyze.getNumRows()), numEigenFaces(0), numFacesAdded(0), numFacesRemoved(0), numOfTopEigenFaces(0) {

    ComplexNum percentVarianceCNum(percentVariance, 0);


    // Make a local copy for processing
    Matrix faceMatrix = faceMatrixToAnalyze;
    // Proceed with scaling and other operations
    faceMatrix = faceMatrix * (1.0 / 255);

    // Assign the scaled matrix to faceDatabase
    faceDatabase = faceMatrix;


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

    while ((eigenSum/totalVariance).getRealPart() < (percentVarianceCNum).getRealPart() && k < eigenvalues.size()) {
        eigenSum += eigenvalues[k];
        k += 1;
    }

    // Build face matrix of eigenfaces
    Matrix eigenFaceDatabaseBuilder(eigenFaces.size(), meanMatrixTranspose.getNumRows()); // Double-check dimensions are correct



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


Matrix FaceSpace::matchFace(Matrix &faceToMatch, double percentSensitivity, int& faceIndex) {

    if (percentSensitivity > 1 || percentSensitivity < 0) {
        throw std::invalid_argument("FaceSpace::matchFace: percentSensitivity must be between 0 and 1 (inclusive)!");
    }

    Matrix scaledFaceToMatch = faceToMatch * (1.0/255);

    Matrix faceDatabaseCopy = faceDatabase;
    int originalFaceNumRows = scaledFaceToMatch.getNumRows();
    int originalFaceNumCols = scaledFaceToMatch.getNumCols();
    Matrix meanFaceMatrix = meanFace;


    Matrix faceToMatchFlattened = flattenMatrix(scaledFaceToMatch);

    Matrix faceToMatchColumnVector = conjTranspose(&faceToMatchFlattened);
    Matrix columnMeanFace = conjTranspose(&meanFaceMatrix);

    std::cout << "faceToMatchColumnVector dimensions: " << faceToMatchColumnVector.getNumRows() << " x " << faceToMatchColumnVector.getNumCols() << std::endl;
    std::cout << "columnMeanFace dimensions: " << columnMeanFace.getNumRows() << " x " << columnMeanFace.getNumCols() << std::endl;


    Matrix processedFaceColumnVector = faceToMatchColumnVector - columnMeanFace;

    Matrix faceToMatchWeightVector(numOfTopEigenFaces, 1);


    for (int i = 0; i < numOfTopEigenFaces; i++) {
        Matrix flattenedEigenFace(1, eigenFaceDatabase.getNumCols());

        for (int j = 0; j < eigenFaceDatabase.getNumCols(); j++) {
            flattenedEigenFace(0, j) = eigenFaceDatabase(i, j);
        }


        Matrix columnEigenFace = conjTranspose(&flattenedEigenFace);
        ComplexNum projectionResult = innerProduct(columnEigenFace, processedFaceColumnVector);
        faceToMatchWeightVector(i, 0) = projectionResult;
    }

    normalizeVectorsInMatrix(&faceToMatchWeightVector);


    // Now go through each vector, project it to the eigenspace, and compare

    std::cout << "going through each vector:" << std::endl;

    double minDistance = 1.5;
    int minDistFaceIndex = 0;



    for (int i = 0; i < numFaces; i++) {

        Matrix flattenedFace(1, faceDatabase.getNumCols());

        for (int k = 0; k < faceDatabase.getNumCols(); k++) {
            flattenedFace(0, k) = faceDatabase(i, k);
        }

        flattenedFace = flattenedFace - meanFace;

        Matrix projectionVector(numOfTopEigenFaces, 1);

        for (int j = 0; j < numOfTopEigenFaces; j++) {

            Matrix columnFace = conjTranspose(&flattenedFace);
            Matrix flattenedEigenFace(1, eigenFaceDatabase.getNumCols());

            for (int k = 0; k < eigenFaceDatabase.getNumCols(); k++) {
                flattenedEigenFace(0, k) = eigenFaceDatabase(j, k);
            }

            Matrix columnEigenFace = conjTranspose(&flattenedEigenFace);
            ComplexNum projectionResult = innerProduct(columnEigenFace, columnFace);
            projectionVector(j, 0) = projectionResult;
        }

        normalizeVectorsInMatrix(&projectionVector);





        Matrix difference = projectionVector - faceToMatchWeightVector;

        double distance = sqrt(VectorNorm(&difference));
        std::cout << "Candidate distance for face " << i << " is " << distance << std::endl;


        if (distance < minDistance) {
            minDistance = distance;
            minDistFaceIndex = i;
        }

    }

    double thresholdDistance = (1 - percentSensitivity) * sqrt(2);


    if (percentSensitivity == 1.0 && minDistance == 0) {
        std::cout << "Perfect match found!" << std::endl;
        faceIndex = minDistFaceIndex;
        Matrix match(1, faceDatabase.getNumCols());

        for (int k = 0; k < faceDatabase.getNumCols(); k++) {
            ComplexNum cNumb = faceDatabaseCopy(minDistFaceIndex, k);
            match(0, k) = cNumb;
        }

        Matrix reconstructedFace = unflattenMatrix(match, originalFaceNumRows, originalFaceNumCols);
        return reconstructedFace;
    }

    else if (minDistance > thresholdDistance) {
        std::cout << "FaceSpace::matchFace: No match found!" << std::endl;
        faceIndex = 0;
        return {1, 1};
    } else {
        std::cout << "Match found with face index: " << minDistFaceIndex << std::endl;
        Matrix match(1, faceDatabase.getNumCols());

        for (int k = 0; k < faceDatabase.getNumCols(); k++) {
            ComplexNum cNumb = faceDatabaseCopy(minDistFaceIndex, k);
            match(0, k) = cNumb;
        }

        faceIndex = minDistFaceIndex;
        Matrix reconstructedFace = unflattenMatrix(match, originalFaceNumRows, originalFaceNumCols);
        return reconstructedFace;
    }



}





