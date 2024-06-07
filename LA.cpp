#include "LA.h"


Matrix matMul(Matrix* lhs, Matrix* rhs) {
    if (lhs->getNumCols() != rhs->getNumRows()) {
        throw std::invalid_argument("Invalid matrix dimensions");
    }

    Matrix product(lhs->getNumRows(), rhs->getNumCols());
    for (int i = 0; i < lhs->getNumRows(); i++) {
        for (int j = 0; j < rhs->getNumCols(); j++) {
            ComplexNum sum(0, 0);
            for (int k = 0; k < lhs->getNumCols(); k++) {
                sum = sum + ((*lhs)(i,k) * (*rhs)(k,j));
            }
            product(i,j) = sum;
        }
    }

    return product;
}

double realTraceOfMatrix(Matrix* matrixToTrace) { // Only defined for real matrices at the moment
    if (matrixToTrace->getNumRows() != matrixToTrace->getNumCols()) {
        throw std::invalid_argument("Trace is only defined for square matrices");
    }

    double runningTrace = 0;
    for (int i = 0; i < matrixToTrace->getNumRows(); i++) {
        double numberToAddToRunningTrace = complexNumToDouble((*matrixToTrace)(i,i));
        runningTrace += numberToAddToRunningTrace;
    }

    return runningTrace;
}

Matrix conjTranspose(Matrix* matrixToTranspose) {
    Matrix matrixToReturn(matrixToTranspose->getNumCols(), matrixToTranspose->getNumRows());
    for (int i = 0; i < matrixToTranspose->getNumRows(); i++) {
        for (int j = 0; j < matrixToTranspose->getNumCols(); j++) {
            ComplexNum transposedEntry = complexConjugate((*matrixToTranspose)(i,j));
            matrixToReturn(j,i) = transposedEntry;
        }
    }
    return matrixToReturn;
}

double frobeniusNorm(Matrix* matrixToNorm) {
    Matrix conjugateTransposeMatrix = conjTranspose(matrixToNorm);
    Matrix matrixStarMatrix = matMul(&conjugateTransposeMatrix, matrixToNorm);
    //std::cout << "Matrix mul has dimensions " << matrixStarMatrix.numRows << "x" << matrixStarMatrix.numCols << std::endl;
    double trace = realTraceOfMatrix(&matrixStarMatrix);
    return sqrt(trace);
}

ComplexNum InnerProduct(Matrix& u, Matrix& v) {
    if (u.getNumCols() != 1 || v.getNumCols() != 1) {
        throw std::invalid_argument("Vectors must be nx1 matrices");
    }

    ComplexNum result(0, 0);
    for (int i = 0; i < u.getNumRows(); i++) {
        ComplexNum vEntryConjugate = complexConjugate(v(i, 0));
        result += (u(i, 0) * vEntryConjugate);
    }
    return result;
}


// Method to normalize vectors in a matrix (alters matrix upon which it is called)
void normalizeVectorsInMatrix(Matrix* pointerToMatrix) {
    for (int i = 0; i < pointerToMatrix->getNumCols(); i++) {

        double Norm = 0;
        for (int j = 0; j < pointerToMatrix->getNumRows(); j++) {
            double magnitude = magnitudeOfNumber((*pointerToMatrix)(j,i));
            Norm += (magnitude * magnitude);
        }

        Norm = sqrt(Norm);
        double oneOverMagnitude = 1/Norm;
        for (int j = 0; j < pointerToMatrix->getNumRows(); j++) {
            (*pointerToMatrix)(j,i) = (*pointerToMatrix)(j,i) * oneOverMagnitude;
        }

    }
}

Matrix GramSchmidt(Matrix const& M) {
    Matrix result(M.getNumRows(),M.getNumCols());

    for(int i = 0; i < M.getNumCols(); i++) {
        Matrix Vk = M[i];
        Matrix Uk = Vk;
        for(int j = 0; j < i; j++) {

            Matrix Un = result[j];
            ComplexNum VdotU = InnerProduct(Vk,Un);
            ComplexNum UdotU = InnerProduct(Un,Un);

            //compute the projection
            ComplexNum quotient = (VdotU / UdotU);
            Un = Un * quotient;
            Uk = Uk - Un;
        }
        result.columnAssign(i,&Uk);
    }

    normalizeVectorsInMatrix(&result);
    return result;
}

void swapRowsInMatrix(Matrix* matrixToSwap, int posRowOne, int posRowTwo) {

    Matrix temp(1, matrixToSwap->getNumCols());
    for (int i = 0; i < matrixToSwap->getNumCols(); i++) {
        temp(0, i) = (*matrixToSwap)(posRowOne, i);
    }

    for (int i = 0; i < matrixToSwap->getNumCols(); i++) {
        (*matrixToSwap)(posRowOne, i) = (*matrixToSwap)(posRowTwo, i);
    }

    for (int i = 0; i < matrixToSwap->getNumCols(); i++) {
        (*matrixToSwap)(posRowTwo, i) = temp(0, i);
    }
}

Matrix createAugmentedMatrix(Matrix* matrixA, Matrix* matrixB) {
    if (matrixA->getNumRows() != matrixB->getNumRows()) {
        throw std::invalid_argument("Matrices must have same number of rows");
    }
    Matrix matrixToReturn(matrixA->getNumRows(), matrixA->getNumCols() + matrixB->getNumCols());

    for (int i = 0; i < matrixA->getNumRows(); i++) {
        for (int j = 0; j < matrixA->getNumCols(); j++) {
            matrixToReturn(i, j) = (*matrixA)(i, j);
        }
        for (int j = 0; j < matrixB->getNumCols(); j++) {
            matrixToReturn(i, j + matrixA->getNumCols()) = (*matrixB)(i, j);
        }
    }
    return matrixToReturn;
}

// Row reduction to check for linear independence will be a separate method called "row reduction"
// Will return a vector later
Matrix gaussianElimination(Matrix* matrix, Matrix* vector) {

    Matrix augmentedMatrixPrePointer = createAugmentedMatrix(matrix, vector); // Temporary fix
    Matrix* augmentedMatrix = &augmentedMatrixPrePointer;

    if (augmentedMatrix->getNumRows() + 1 != augmentedMatrix->getNumCols()) {
        throw std::invalid_argument("Augmented matrix must be of size n x n+1");
    }
    // Forward elimination
    for (int i = 0; i < augmentedMatrix->getNumCols() - 2; i++) {
        // Get the largest element in column and swap rows for numerical stability

        ComplexNum maxElementInColumn(0,0);
        int maxElementInColumnRowPos = i;
        for (int j = i; j < augmentedMatrix->getNumRows(); j++) {
            if (magnitudeOfNumber((*augmentedMatrix)(j, i)) > magnitudeOfNumber(maxElementInColumn)) {
                maxElementInColumnRowPos = j;
                maxElementInColumn = (*augmentedMatrix)(j, i);
            }
        }
        if (magnitudeOfNumber(maxElementInColumn) < 1e-6) {
            throw std::invalid_argument("No non-zero element in column");
        }

        if (maxElementInColumnRowPos != i) {
            swapRowsInMatrix(augmentedMatrix, i, maxElementInColumnRowPos);
        }

        //std::cout << "LARGEST ELEMENT IN COLUMN IS " << maxElementInColumn << std::endl;

        // Use operations to ensure the column is consistent with a matrix in REF

        for (int j = i + 1; j < augmentedMatrix->getNumRows(); j++) {
            // Add multiples of j = i to each row from i + 1 to the end of column to get zeros
            ComplexNum scalarMultipleForRowi = ((*augmentedMatrix)(j, i))/(*augmentedMatrix)(i, i);
            scalarMultipleForRowi = scalarMultipleForRowi * -1;

            //std::cout << "scalar multiple for row " << i + 1 << " is " << scalarMultipleForRowi << std::endl;

            for (int k = 0; k < augmentedMatrix->getNumCols(); k++) {
                (*augmentedMatrix)(j, k) = (*augmentedMatrix)(j, k) + ((*augmentedMatrix)(j - 1, k) * scalarMultipleForRowi);
            }

        }

    }

    // Execute the same procedure now considering the entries above the diagonal (this can be replaced with straight back substitution in the future)
    for (int i = 0; i < augmentedMatrix->getNumCols() - 1; i++) {
        // Get the largest element in column and swap rows for numerical stability
        if (magnitudeOfNumber((*augmentedMatrix)(i, i)) < 1e-6) {
            throw std::invalid_argument("Matrix is singular!");
        }
        //std::cout << "here?" << std::endl;

        for (int j = i - 1; j >= 0; j--) {

            // Add multiples of j = i to each row from i + 1 to the end of column to get zeros
            ComplexNum scalarMultipleForRowi = ((*augmentedMatrix)(j, i))/(*augmentedMatrix)(i, i);
            scalarMultipleForRowi = scalarMultipleForRowi * -1;

            //std::cout << "scalar multiple for row " << i + 1 << " is " << scalarMultipleForRowi << std::endl;

            for (int k = 0; k < augmentedMatrix->getNumCols(); k++) {
                (*augmentedMatrix)(j, k) = (*augmentedMatrix)(j, k) + ((*augmentedMatrix)(j + 1, k) * scalarMultipleForRowi);
            }

        }

    }

    //std::cout << augmentedMatrixPrePointer << std::endl;
    for (int j = 0; j < augmentedMatrix->getNumRows(); j++) {
        ComplexNum reciprocal = ComplexNum(1, 0) / (*augmentedMatrix)(j,j);
        for (int k = 0; k < augmentedMatrix->getNumCols(); k++) {
            (*augmentedMatrix)(j,k) = (*augmentedMatrix)(j,k) * reciprocal;
        }
    }

    Matrix vectorToReturn(vector->getNumRows(), vector->getNumCols());
    for (int j = 0; j < vector->getNumRows(); j++) {
        vectorToReturn(j, 0) = (*augmentedMatrix)(j, augmentedMatrix->getNumCols() - 1);
    }

    return vectorToReturn;
}

Matrix inverseMatrix(Matrix* matrixToInvert) {
    if (matrixToInvert->getNumRows() != matrixToInvert->getNumCols()) {
        throw std::invalid_argument("Matrix must be n x n!");
    }

    Matrix outputMatrix(matrixToInvert->getNumRows(), matrixToInvert->getNumCols());
    for (int i = 0; i < matrixToInvert->getNumRows(); i++) {
        Matrix unitVector(matrixToInvert->getNumRows(), 1);
        unitVector(i, 0) = 1;
        Matrix resultVector = gaussianElimination(matrixToInvert, &unitVector);

        for (int j = 0; j < matrixToInvert->getNumRows(); j++) {
            outputMatrix(j, i) = resultVector(j , 0);
        }
    }
    return outputMatrix;
}




