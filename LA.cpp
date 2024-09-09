#include "LA.h"
#include "Matrix.h"
#include <opencv2/opencv.hpp>
#include "data.h"

void thresholdStabilize(Matrix* matrixToStabilize) {
    for (int i = 0; i < matrixToStabilize->getNumRows(); i++) {
        for (int j = 0; j < matrixToStabilize->getNumCols(); j++) {
            if (abs((*matrixToStabilize)(i,j).getRealPart()) < 1e-12) {
                (*matrixToStabilize)(i,j) = ComplexNum(0, (*matrixToStabilize)(i,j).getImagPart());
                //std::cout << "stab called";
            }
            if (abs((*matrixToStabilize)(i,j).getImagPart()) < 1e-12) {
                (*matrixToStabilize)(i,j) = ComplexNum((*matrixToStabilize)(i,j).getRealPart(), 0);
                //std::cout<< "stab called";
            }

        }
    }
}


Matrix matMul(Matrix* lhs, Matrix* rhs) {
    if (lhs->getNumCols() != rhs->getNumRows()) {
        throw std::invalid_argument("matMul: Invalid matrix dimensions");
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
    thresholdStabilize(&product);
    return product;
}

//element wise multiplication
Matrix hadamardProduct(Matrix* lhs, Matrix* rhs) {
    if(lhs->getNumCols() != rhs->getNumCols() || lhs->getNumRows() != rhs->getNumRows()) {
        throw std::invalid_argument("hadamardProduct: Matrices have incongruent dimensions");
    }

    Matrix result = *lhs;
    for(int i = 0; i < lhs->getNumRows(); i++) {
        for(int j = 0; j < lhs->getNumCols(); j++) {
            result(i,j) = (*lhs)(i,j) * (*rhs)(i,j);
        }
    }
    return result;
}

double realTraceOfMatrix(Matrix* matrixToTrace) { // Only defined for real matrices at the moment
    if (matrixToTrace->getNumRows() != matrixToTrace->getNumCols()) {
        throw std::invalid_argument("realTraceOfMatrix: Trace is only defined for square matrices");
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

Matrix transpose(Matrix* matrixToTranspose) {
    Matrix matrixToReturn(matrixToTranspose->getNumCols(), matrixToTranspose->getNumRows());
    for (int i = 0; i < matrixToTranspose->getNumRows(); i++) {
        for (int j = 0; j < matrixToTranspose->getNumCols(); j++) {
            ComplexNum transposedEntry = (*matrixToTranspose)(i,j);
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

ComplexNum innerProduct(Matrix& u, Matrix& v) {
    if (u.getNumCols() != 1 || v.getNumCols() != 1) {
        throw std::invalid_argument("InnerProduct: Vectors must be nx1 matrices");
    }

    ComplexNum result(0, 0);
    for (int i = 0; i < u.getNumRows(); i++) {
        ComplexNum vEntryConjugate = complexConjugate(v(i, 0));
        result += (u(i, 0) * vEntryConjugate);
    }
    return result;
}

ComplexNum mean(Matrix& u) {
    ComplexNum mean(0,0);
    for(int i = 0; i < u.getNumRows(); i++) {
        for(int j = 0; j < u.getNumCols(); j++) {
            mean += u(i,j);
        }
    }
    return (mean / (u.getNumCols() * u.getNumRows()));
}


// Method to normalize vectors in a matrix (alters matrix upon which it is called)
void normalizeVectorsInMatrix(Matrix* pointerToMatrix) {
    for (int i = 0; i < pointerToMatrix->getNumCols(); i++) {

        double Norm = 0;
        for (int j = 0; j < pointerToMatrix->getNumRows(); j++) {
            double magnitude = magnitudeOfNumber((*pointerToMatrix)(j,i));
            Norm += (magnitude * magnitude);
        }
        if (Norm == 1) {
            return;
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
            ComplexNum VdotU = innerProduct(Vk,Un);
            ComplexNum UdotU = innerProduct(Un,Un);

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
        throw std::invalid_argument("createAugmentedMatrix: Matrices must have same number of rows");
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
Matrix gaussianElimination(Matrix* matrix, Matrix* vector) {

    Matrix augmentedMatrixPrePointer = createAugmentedMatrix(matrix, vector); // Temporary fix
    Matrix* augmentedMatrix = &augmentedMatrixPrePointer;

    //std::cout << "Augmented matrix to perform row reduction on:\n" << augmentedMatrixPrePointer << std::endl;

    if (augmentedMatrix->getNumRows() + 1 != augmentedMatrix->getNumCols()) {
        throw std::invalid_argument("gaussianElimination: Augmented matrix must be of size n x n+1");
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
            throw std::invalid_argument("gaussianElimination: Matrix is singular!");
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

            //std::cout << augmentedMatrixPrePointer << std::endl;

            //std::cout << "scalar multiple for row " << j << " is " << scalarMultipleForRowi << std::endl;

            for (int k = 0; k < augmentedMatrix->getNumCols(); k++) {
                (*augmentedMatrix)(j, k) = (*augmentedMatrix)(j, k) + ((*augmentedMatrix)(i, k) * scalarMultipleForRowi);
                //std::cout << (*augmentedMatrix)(j, k) << std::endl;
            }
            //std::cout << augmentedMatrixPrePointer << std::endl;

        }

    }
    //std::cout << "Upper triangular matrix is\n" << augmentedMatrixPrePointer << "\n" << std::endl;

    // Execute the same procedure now considering the entries above the diagonal (this can be replaced with straight back substitution in the future)
    for (int i = 0; i < augmentedMatrix->getNumCols() - 1; i++) {
        // Get the largest element in column and swap rows for numerical stability
        if (magnitudeOfNumber((*augmentedMatrix)(i, i)) < 1e-6) {
            throw std::invalid_argument("gaussianElimination: Matrix is singular!");
        }
        //std::cout << "here?" << std::endl;

        for (int j = i - 1; j >= 0; j--) {

            // Add multiples of j = i to each row from i + 1 to the end of column to get zeros
            ComplexNum scalarMultipleForRowi = ((*augmentedMatrix)(j, i))/(*augmentedMatrix)(i, i);
            scalarMultipleForRowi = scalarMultipleForRowi * -1;

            //std::cout << "scalar multiple for row " << i + 1 << " is " << scalarMultipleForRowi << std::endl;

            for (int k = 0; k < augmentedMatrix->getNumCols(); k++) {
                (*augmentedMatrix)(j, k) = (*augmentedMatrix)(j, k) + ((*augmentedMatrix)(i, k) * scalarMultipleForRowi);
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

Matrix inverseMatrix(Matrix* matrixToInvert) { // Really, really slow. That being said, we can speed this up in the future if needed, but at the moment I want "easy" code for debugging
    if (matrixToInvert->getNumRows() != matrixToInvert->getNumCols()) {
        throw std::invalid_argument("inverseMatrix: Matrix must be nxn!");
    }

    Matrix outputMatrix(matrixToInvert->getNumRows(), matrixToInvert->getNumCols());
    for (int i = 0; i < matrixToInvert->getNumRows(); i++) {
        Matrix unitVector(matrixToInvert->getNumRows(), 1);
        unitVector(i, 0) = 1;
        Matrix resultVector = gaussianElimination(matrixToInvert, &unitVector);

        for (int j = 0; j < matrixToInvert->getNumRows(); j++) {
            outputMatrix(j, i) = resultVector(j , 0);
        }
        //std::cout << outputMatrix[i] << std::endl;
    }

    return outputMatrix;
}

//gets magnitude of column matrixcolumnFace
double VectorNorm(Matrix* pointerToMatrix) {
    double Norm = 0;
    for (int j = 0; j < pointerToMatrix->getNumRows(); j++) {
        double magnitude = magnitudeOfNumber((*pointerToMatrix)(j,0));
        Norm += (magnitude * magnitude);
    }
    if(Norm == 1) {
        return Norm;
    }
    return sqrt(Norm);
}

Matrix unitVector(int k, int dim) {
    if (k > dim || dim < 1  || k < 0) {
        throw std::invalid_argument("unitVector: Invalid input");
    }
    Matrix e(dim,1);
    e(k,0) = 1;
    return e;
}

Matrix identityMatrix(int dim) {
    if(dim < 1) {
        throw std::invalid_argument("identityMatrix: Invalid input");
    }
    Matrix I(dim,dim);

    //replace each column with unit vector
    for(int i = 0; i < dim; i++) {
        Matrix ek = unitVector(i,dim);
        I.columnAssign(i,&ek);
    }

    return I;
}

Matrix frontFillVec(Matrix M, int dim, ComplexNum const& Fill) {
    if(M.getNumCols() > 1 || dim < 1 || dim < M.getNumCols()) {
        throw std::invalid_argument("frontFillVec: Invalid input");
    }


    Matrix result(dim,1);
    for(int i = 0; i < dim; i++) {
        //fill with num until getting to matrix with data
        if(i < dim - M.getNumRows()) {
            result(i,0) = Fill;
        }
        else {
            result(i,0) = M(i - (dim - M.getNumRows()),0);
        }
    }
    return result;
}

Matrix minorMatrix(Matrix M,int iDel, int jDel) {
    Matrix Minor(M.getNumRows() -1,M.getNumCols() -1);
    for(int i = 0; i < M.getNumRows(); i++) {
        if (i < iDel) {
            for(int j = 0; j < M.getNumCols(); j++) {{
                if(j < jDel) {
                    Minor(i,j) = M(i,j);
                }
                else if(j > jDel) {
                    Minor(i,j-1) = M(i,j);
                }
            }}
        }
        else if(i > iDel) {
            for(int j = 0; j < M.getNumCols(); j++) {{
                if(j < jDel) {
                    Minor(i-1,j) = M(i,j);
                }
                else if(j > jDel) {
                    Minor(i-1,j-1) = M(i,j);
                }
            }}
        }
    }
    return Minor;
}

Matrix tridiagonalizeMatrix(Matrix& matrixToTri) {

    if (matrixToTri.getNumCols() != matrixToTri.getNumRows()) {
        throw std::invalid_argument("tridiagonalizeMatrix: Matrix is not square!");
    }

    Matrix A = matrixToTri;

    for (int i = 0; i < matrixToTri.getNumRows() - 2; i++) { // Loop

        // Step 1: Calculate alpha
            // Get sign
        int alphaSign = getRealSign(A(i+1, i));
            // Get sum
        ComplexNum alphaSum(0, 0);
        for (int j = i + 1; j < A.getNumRows(); j++) {

            ComplexNum a = A(j, i);

            // Get power
            ComplexNum numSquared = a * a;

            alphaSum += numSquared;
        }
        // Square root sum
        double sqrtSum = sqrt(alphaSum.getRealPart());
        double alpha = -1 * alphaSign * sqrtSum;


        // Step 2: Find r
        ComplexNum entry = A(i+1, i);
        double alphaEntry = (entry * alpha).getRealPart();
        double difference = (1.0/2.0) * ((alpha * alpha) - alphaEntry);

        double r = sqrt(difference);

        // Construct vector
        Matrix vector(A.getNumRows(), 1);

        // Set first i positions to 0
        for (int j = 0; j <= i; j++) {
            vector(j, 0) = ComplexNum(0, 0);
        }

        vector(i+1, 0) = (entry - alpha) * (ComplexNum(1, 0)/(2*r));


        for (int j = i+2; j < A.getNumRows(); j++) {
            vector(j, 0) = (A(j, i) * (ComplexNum(1, 0)/(2*r)));
        }

        // Compute P
        Matrix identity = identityMatrix(A.getNumRows());

        Matrix vectorTranspose = conjTranspose(&vector);
        Matrix vectorTransform = matMul(&vector, &vectorTranspose);
        vectorTransform = vectorTransform * 2;

        Matrix P = identity - vectorTransform;

        A = matMul(&A, &P);
        A = matMul(&P, &A);
    }


    return A;
}





Matrix householderReflection(Matrix* x) { //todo stabilize


    if (x->getNumCols() > 1) {
        throw std::invalid_argument("householderTransform: Invalid input");
    }
    Matrix y = (*x)(0,0).sign() * VectorNorm(x) * unitVector(0,x->getNumRows());
    y = *x + y;
    //std::cout << "y dims: " << y.getNumRows() << " x " << y.getNumCols() << std::endl;
    return y;



}

std::vector<Matrix> QRDecomp(Matrix const& A) { //todo stabilize

    Matrix R = A;
    std::vector<Matrix> HTransforms;
    for(int i = 0; i < A.getNumCols(); i++) {
        if(i > 0) {
            R = minorMatrix(R,0,0);
        }


        //get transformation on column
        Matrix a = R[0];
        Matrix v = householderReflection(&a);
        Matrix vT = conjTranspose(&v);


        //get the transformation on the matrix
        Matrix H = identityMatrix(R.getNumRows());
        H = H - ((ComplexNum(2,0) / innerProduct(v,v) * (matMul(&v,&vT)) ));
        HTransforms.push_back(H);

        R = matMul(&H,&R);


    }

    for(int i = 1; i < HTransforms.size(); i++) {
        Matrix resized = identityMatrix(HTransforms[0].getNumCols());
        for(int j = HTransforms[0].getNumCols()-HTransforms[i].getNumCols(); j < HTransforms[0].getNumCols(); j++){
            Matrix toFrontLoad = (HTransforms[i])[j-(HTransforms[0].getNumCols()-HTransforms[i].getNumCols())];
            Matrix processed = frontFillVec(toFrontLoad,HTransforms[0].getNumCols(),ComplexNum(0,0));
            resized.columnAssign(j,&processed);
        }
        HTransforms[i] = resized;
    }

    Matrix Q = HTransforms.back();
    for(int i = HTransforms.size() - 2; i >= 0; i--) {
        Q = matMul(&HTransforms[i],&Q);
    }

    R = A;
    for(int i = 0; i < HTransforms.size(); i++) {
        R = matMul(&HTransforms[i],&R);
    }
    return {Q,R};
}


bool isUpperTriangular(Matrix* mat) {
    for (int i = 0; i < mat->getNumCols(); ++i) {
        for (int j = i + 1; j < mat->getNumRows(); ++j) {
            if (magnitudeOfNumber((*mat)(j, i)) > 1e-6) {
                return false;
            }
        }
    }
    return true;
}

ComplexNum mu(Matrix* matrix) {
    ComplexNum delta = ((*matrix)(0, 0) - (*matrix)(1, 1))/2;
    ComplexNum radicalMult = sqrt((delta * delta) + ((*matrix)(0, 1)) * ((*matrix)(1, 0)));

    double sign;
    if (delta.getRealPart() >= 0) {
        sign = 1;
    } else {
        sign = -1;
    }

    ComplexNum WilkinsonShift = delta - ((sign) * radicalMult);

    return WilkinsonShift;
}

std::vector<ComplexNum> eigenvalues(Matrix* matrix) {
    if (matrix->getNumRows() != matrix->getNumCols()) {
        throw std::invalid_argument("eigenvalues: matrix must be n x n!");
    }

    if (matrix->getNumCols() == 1 && matrix->getNumRows() == 1) { // Trivial case
        std::vector<ComplexNum> ret;
        ret.push_back(ComplexNum((*matrix)(0, 0)));
        return ret;
    }

    Matrix matrixToIterate = tridiagonalizeMatrix(*matrix);
    // Check to see if matrix is already upper triangular


    bool akIsUpperTriangular = isUpperTriangular(&matrixToIterate);


    while (!akIsUpperTriangular) {

        Matrix subMatrix(2,2);
        int numRows = matrix->getNumRows();
        subMatrix(0, 0) = (*matrix)(numRows - 2, numRows - 2);
        subMatrix(1, 0) = (*matrix)(numRows - 1, numRows - 2);
        subMatrix(0, 1) = (*matrix)(numRows - 2, numRows - 1);
        subMatrix(1, 1) = (*matrix)(numRows - 1, numRows - 1);

        ComplexNum shift = mu(&subMatrix);

        Matrix identity = identityMatrix(matrix->getNumRows());
        matrixToIterate = matrixToIterate - (shift * identity);


        //std::cout << "Matrix to iterate:\n" << matrixToIterate << std::endl;
        std::vector<Matrix> currentQR = QRDecomp(matrixToIterate);
        //std::cout << "Q is:\n" << currentQR[0] << std::endl;
        //std::cout << "R is:\n" << currentQR[1] << std::endl;
        matrixToIterate = matMul(&currentQR[1], &currentQR[0]);
        //std::cout << matrixToIterate << std::endl;


        matrixToIterate = matrixToIterate + (shift * identity);

        // Check if matrixToIterate is upper triangular
        akIsUpperTriangular = isUpperTriangular(&matrixToIterate);
    }
    //std::cout << "akIsUpperTriangular and ak is:\n" << matrixToIterate << std::endl;

    std::set<ComplexNum> setOfEigenvalues;
    for (int i = 0; i < matrixToIterate.getNumRows(); i++) {
        double numImagPart = (matrixToIterate(i, i)).getImagPart();
        double numRealPart = (matrixToIterate(i, i)).getRealPart();
        ComplexNum newCNum(numRealPart, numImagPart);
        setOfEigenvalues.insert(newCNum);
    }

    std::vector<ComplexNum> vectorOfEigenvalues;
    for (const ComplexNum& element : setOfEigenvalues) {
        double numImagPart = element.getImagPart();
        double numRealPart = element.getRealPart();
        ComplexNum newCNum(numRealPart, numImagPart);
        vectorOfEigenvalues.push_back(newCNum);
    }

    return vectorOfEigenvalues;

}


std::vector<Matrix> eigenvectors(Matrix* matrix, std::vector<ComplexNum>& correspondingEigenValues) { // Only works for real eigenvalues

    if (matrix->getNumRows() != matrix->getNumCols()) {
        throw std::invalid_argument("eigenvalues: matrix must be n x n!");
    }

    if (matrix->getNumCols() == 1 && matrix->getNumRows() == 1) { // Trivial case
        std::vector<Matrix> ret;
        Matrix matret(1, 1);
        matret(0, 0) = 1;
        ret.push_back(matret);
        return ret;
    }

    Matrix matrixToIterate = tridiagonalizeMatrix(*matrix);
    // Check to see if matrix is already upper triangular

    Matrix cumulativeQ = identityMatrix(matrix->getNumCols());

    bool akIsUpperTriangular = isUpperTriangular(&matrixToIterate);


    int counter = 0;
    while (!akIsUpperTriangular) {
        counter += 1;

        Matrix subMatrix(2,2);
        int numRows = matrix->getNumRows();
        subMatrix(0, 0) = (*matrix)(numRows - 2, numRows - 2);
        subMatrix(1, 0) = (*matrix)(numRows - 1, numRows - 2);
        subMatrix(0, 1) = (*matrix)(numRows - 2, numRows - 1);
        subMatrix(1, 1) = (*matrix)(numRows - 1, numRows - 1);

        ComplexNum shift = mu(&subMatrix);

        Matrix identity = identityMatrix(matrix->getNumRows());
        matrixToIterate = matrixToIterate - (shift * identity);


        std::vector<Matrix> currentQR = QRDecomp(matrixToIterate);
        matrixToIterate = matMul(&currentQR[1], &currentQR[0]);

        cumulativeQ = matMul(&cumulativeQ, &currentQR[0]);

        matrixToIterate = matrixToIterate + (shift * identity);

        // Check if matrixToIterate is upper triangular
        cumulativeQ = GramSchmidt(cumulativeQ);
        akIsUpperTriangular = isUpperTriangular(&matrixToIterate);
    }
    //std::cout << "akIsUpperTriangular and ak is:\n" << matrixToIterate << std::endl;

    std::vector<Matrix> vectorOfEigenvectors; // Ordered numerically
    std::vector<ComplexNum> vectorOfEigenvalues;
    for (int i = 0; i < matrixToIterate.getNumCols(); i++) {
        vectorOfEigenvalues.push_back( matrixToIterate(i, i));
        vectorOfEigenvectors.push_back(cumulativeQ[i]);
    }



    // Now sort by eigenvalue
    //std::cout << "cum q:\n" << cumulativeQ << std::endl;


    std::vector<Matrix> orderedEigenvectors;
    std::vector<ComplexNum> orderedEigenvalues;
    int vecLength = vectorOfEigenvalues.size();


    for (int i = 0; i < vecLength; i++) {

        double max = magnitudeOfNumber(vectorOfEigenvalues[0]);
        ComplexNum maxEigenVal = vectorOfEigenvalues[0];
        int maxPos = 0;
        for (int j = 0; j < vectorOfEigenvalues.size(); j++) {
            //std::cout << "Vector of eigenvalues size is: " << vectorOfEigenvalues.size() << std::endl;

            if (magnitudeOfNumber(vectorOfEigenvalues[j]) > max) {
                max = magnitudeOfNumber(vectorOfEigenvalues[j]);
                maxPos = j;
                maxEigenVal = vectorOfEigenvalues[j];
            }
        }

        orderedEigenvalues.push_back(maxEigenVal);
        orderedEigenvectors.push_back(vectorOfEigenvectors[maxPos]);

        vectorOfEigenvalues.erase(vectorOfEigenvalues.begin() + maxPos);
        vectorOfEigenvectors.erase(vectorOfEigenvectors.begin() + maxPos);
    }

    // We now have two matching vectors ordered from largest eigenvalue to smallest.

    /*
    for (int i = 0; i < orderedEigenvalues.size(); i++) {
        std::cout << orderedEigenvalues[i] << std::endl;
        std::cout << orderedEigenvectors[i] << std::endl;
    }
    */

    correspondingEigenValues = orderedEigenvalues;
    return orderedEigenvectors;



}









std::vector<Matrix> singularValueDecomp(Matrix* matrix) {
    // The first step is to create an empty sigma matrix, as this part of the code finds the singular values and sigma
    std::vector<Matrix> decompToReturn;

    Matrix sigma(matrix->getNumRows(), matrix->getNumCols());

    // Now create matrix * matrix^T and matrix^T * matrix
    Matrix matrixConjTranspose = conjTranspose(matrix);
    Matrix AAT = matMul(matrix, &matrixConjTranspose);
    Matrix ATA = matMul(&matrixConjTranspose, matrix);

    // One can choose which matrix to use for eigenvalues based on the dimensions of the original matrix if needed, but only do this if speed becomes a legitimate concern
    // We now want to find singular values

    std::vector<ComplexNum> eigenvals = eigenvalues(&AAT);

    std::vector<ComplexNum> singularVals;
    //singularVals.reserve(eigenvals.size());
    for (int i = 0; i < eigenvals.size(); i++) {
        singularVals.push_back(sqrt(eigenvals[i]));
    }


    // Now we construct sigma by placing the singular values in order
    int iCounter = 0;
    for (int i = singularVals.size() - 1; i >= 0; i--) {
        if (iCounter >= matrix->getNumRows() || iCounter >= matrix->getNumCols()) {
            break;
        }
        sigma(iCounter, iCounter) = ComplexNum(singularVals[i]);
        iCounter += 1;
    }

    //std::cout << sigma;

    // We now proceed to construct V

    std::vector<ComplexNum> corresEigensV;
    std::vector<Matrix> vectorsInV = eigenvectors(&ATA, corresEigensV);

    Matrix V(matrix->getNumCols(), matrix->getNumCols());
    for (int i = 0; i < vectorsInV.size(); i++) {
        for (int j = 0; j < V.getNumRows(); j++) {
            V(j, i) = (vectorsInV[i])(j, 0);
        }
    }

    //std::cout << "V:\n" << V << std::endl;


    // Now to construct U
    std::vector<ComplexNum> corresEigensU;
    std::vector<Matrix> vectorsInU = eigenvectors(&AAT, corresEigensU);


    Matrix U(matrix->getNumRows(), matrix->getNumRows());
    for (int i = 0; i < vectorsInU.size(); i++) {
        //U[i] = vectorsInU[i];
        //std::cout << "here?" << std::endl;
        for (int j = 0; j < U.getNumRows(); j++) {
            //kstd::cout << "accessing (" << j << "," << i << ")" << std::endl;
            U(j, i) = (vectorsInU[i])(j, 0);
        }
    }

    // Run algorithm for partial reconstruction and adjust signs


    /*
    for (int i = 0; i < sigma.getNumRows(); i++) {

        Matrix rightSide = VTranpose(i);
        Matrix leftSide = U[i];
        Matrix partialReconstruction = matMul(&rightSide, &leftSide) * sigma(i, i);


        Matrix errorEst1 = (partialReconstruction - (matMul(&rightSide, &leftSide) * sigma(i, i)));
        Matrix errorEst2 = (-1 * partialReconstruction - (matMul(&rightSide, &leftSide) * sigma(i, i)));
        //Matrix errorEst3 =

        double errorOriginal = frobeniusNorm(&errorEst1);
        double errorFlipped = frobeniusNorm(&errorEst2);

        if (errorFlipped < errorOriginal) { // Needs generealization
            U[i] = U[i] * -1;
        }

    }

     */

    // These are rough drafts of various sign-flipping algorithms. Do not delete them; these algorihms might prove useful in the future if the dot product algorithm turns out to not be robust

    /*
    Matrix Sigma_inv(matrix->getNumCols(), matrix->getNumCols());
    for (int i = 0; i < matrix->getNumCols(); ++i) {
        if (sigma(i, i).getRealPart() != 0) { // Assuming ComplexNum has a member 'real'
            Sigma_inv(i, i) = 1.0 / sigma(i, i).getRealPart(); // Assuming ComplexNum supports double division
        }
    }


    Matrix temp = matMul(&V, &Sigma_inv);
    U = matMul(matrix, &temp);
    //U = A * V * Sigma_inv;
     */


    int colMin = V.getNumCols();
    if (U.getNumCols() < colMin) {
        colMin = U.getNumCols();
    }

    for (int i = 0; i < colMin; i++) {
        Matrix col = V[i];
        Matrix Av = matMul(matrix, &col);
        Matrix uCol = U[i];
        ComplexNum dotprodascnum = innerProduct(uCol, Av);
        double dotprod = dotprodascnum.getRealPart();

        if (dotprod < 0) {
             Matrix uColFlipped = -1 * uCol;
             for (int j = 0; j < uColFlipped.getNumRows(); j++) {
                 U(j, i) = uColFlipped(j, 0);
             }
        }
    };



    Matrix VTranpose = transpose(&V);
    decompToReturn.push_back(U);
    decompToReturn.push_back(sigma);
    decompToReturn.push_back(VTranpose);


    return decompToReturn;
}





//assumes all values in the matrix have uniform weight (add another functions that uses weight vector??)
ComplexNum expectedValue(Matrix* Z) {
    if(Z->getNumRows() > 1) {
        throw std::invalid_argument("expectedValue: Not a complex random variable");
    }

    ComplexNum E(0,0);
    for(int i = 0; i < Z->getNumCols(); i++) {
        E += (*Z)(0,i);
    }
    E = E * (1.0/(Z->getNumCols()));
    return E;
}

ComplexNum covariance(Matrix* Z, Matrix* W) {
    Matrix Ztranspose = transpose(Z);
    Matrix Wtranspose = transpose(W);
    Matrix conjW = Wtranspose.conjugate();
    Matrix ZconjW = hadamardProduct(&Ztranspose,&conjW);
    ComplexNum EZconjW = expectedValue(&ZconjW);

    ComplexNum EZ = expectedValue(Z);
    ComplexNum EW = expectedValue(W);

    return EZconjW - (EZ * EW);
}

Matrix covarianceMatrix(Matrix* M) {
    Matrix covarianceMatrix(M->getNumRows(),M->getNumCols());

    //compute mean of each column and subtract to center

    for(int i = 0; i < M->getNumCols(); i++) {
        Matrix column = (*M)[i];
        ComplexNum Mean = mean(column);

        Matrix centered = column - Mean;
        covarianceMatrix.columnAssign(i,&centered);
    }

    Matrix adjoint = conjTranspose(&covarianceMatrix);
   // std::cout<< covarianceMatrix << std::endl << adjoint;
    covarianceMatrix = matMul(&adjoint,&covarianceMatrix) / ComplexNum((M->getNumCols() - 1)*2,0) ;
    return covarianceMatrix;
}


Matrix littleCovariance(Matrix& matrix) {
    //std::cout << matrix << std::endl;
    Matrix matrixTranspose = conjTranspose(&matrix);
    ComplexNum cNum(1.0/(matrix.getNumRows() - 1), 0);
    matrixTranspose = matrixTranspose * cNum;
    Matrix product = matMul(&matrix, &matrixTranspose);
    //std::cout << product << std::endl;
    //product = product * (1.0/(product.getNumRows() - 1));
    return product;
}


std::vector<ComplexNum> CVEigenValues(Matrix& matrix) {
    std::vector<ComplexNum> vectorToReturn;
    cv::Mat CVMatrix = convertMatrixToCVMatrix(matrix);
    cv::Mat CVEigenValues;

    cv::eigen(CVMatrix, CVEigenValues);

    int numEigenValues = CVEigenValues.rows;

    for (int i = 0; i < numEigenValues; i++) {
        double realValue = CVEigenValues.at<double>(i, 0); // Eigenvalue at row i, column 0
        ComplexNum Cnum(realValue, 0);
        vectorToReturn.push_back(Cnum);
    }

    return vectorToReturn;
}

std::vector<Matrix> CVEigenvectors(Matrix* matrix, std::vector<ComplexNum>& correspondingEigenValues) {
    std::vector<ComplexNum> eigenValues;
    std::vector<Matrix> eigenVectorMatrices;

    cv::Mat CVMatrix = convertMatrixToCVMatrix(*matrix);
    cv::Mat CVEigenValues;
    cv::Mat CVEigenVectors;

    cv::eigen(CVMatrix, CVEigenValues, CVEigenVectors);
    int numEigenVectors = CVEigenVectors.rows;

    // Extract eigenvectors
    for (int i = 0; i < numEigenVectors; ++i) {
        Matrix eigenVectorMatrix(numEigenVectors, 1);

        for (int j = 0; j < numEigenVectors; ++j) {
            double value = CVEigenVectors.at<double>(i, j);
            eigenVectorMatrix(j, 0) = ComplexNum(value, 0);
        }

        eigenVectorMatrices.push_back(eigenVectorMatrix);
    }


    // Extract eigenvalues
    for (int i = 0; i < numEigenVectors; i++) {
        double realValue = CVEigenValues.at<double>(i, 0); // Eigenvalue at row i, column 0
        ComplexNum Cnum(realValue, 0);
        eigenValues.push_back(Cnum);
    }

    correspondingEigenValues = eigenValues;
    return eigenVectorMatrices;

}


double smallestNumberInRealMatrix(Matrix& matrix) {
    double currentSmallest = matrix(0, 0).getRealPart();
    for (int i = 0; i < matrix.getNumRows(); i++) {
        for (int j = 0; j < matrix.getNumCols(); j++) {
            if (matrix(i, j).getRealPart() < currentSmallest) {
                currentSmallest = matrix(i, j).getRealPart();
            }
        }
    }

    return currentSmallest;
}

double largestNumberInRealMatrix(Matrix& matrix) {
    double currentLargest = matrix(0, 0).getRealPart();
    for (int i = 0; i < matrix.getNumRows(); i++) {
        for (int j = 0; j < matrix.getNumCols(); j++) {
            if (matrix(i, j).getRealPart() > currentLargest) {
                currentLargest = matrix(i, j).getRealPart();
            }
        }
    }

    return currentLargest;
}

Matrix normalizeMatrix(Matrix& matrix, double maxVal, double minVal) {
    double difference = maxVal - minVal;
    Matrix normedMatrix(matrix.getNumRows(), matrix.getNumCols());

    for (int i = 0; i < matrix.getNumRows(); i++) {
        for (int j = 0; j < matrix.getNumCols(); j++) {
            double normedVal = (matrix(i, j).getRealPart() - minVal) / (difference) * 255.0;

            normedMatrix(i, j) = static_cast<unsigned char>(std::round(std::max(0.0, std::min(255.0, normedVal))));
        }

    }

    return normedMatrix;

}