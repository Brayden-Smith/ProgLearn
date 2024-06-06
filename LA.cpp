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
Matrix gaussianElimination(Matrix* augmentedMatrix, Matrix* vector) {
    if (augmentedMatrix->getNumRows() != augmentedMatrix->getNumCols() + 1) {
        throw std::invalid_argument("Augmented matrix must be of size n x n+1");
    }
    // Forward elimination
    for (int i = 0; i < augmentedMatrix->getNumCols() - 2; i++) {
        // Get largest element in column for numerical stability




    }

    // Ensure diagonal entries are non-zero

    // Back substitution

    //
}




