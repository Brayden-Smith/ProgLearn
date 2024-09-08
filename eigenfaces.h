#include "data.h"

Matrix meanFlattenedFace(Matrix matrix);

class FaceSpace {
private:
    Matrix faceDatabase;
    Matrix meanFace;
    Matrix eigenFaceDatabase;
    long numFaces;
    unsigned int numEigenFaces;
    unsigned int numOfTopEigenFaces;

    // For adding and removing faces (will come later)
    unsigned int numFacesAdded;
    unsigned int numFacesRemoved;

public:
    // Constructors
    FaceSpace(Matrix faceMatrixToAnalyze, double percentVariance); // The percentVariance parameter allows the user to select what percent of the variance of faces the eigenfaces should cover
    FaceSpace(std::string& path, double percentVariance);

    long getNumFaces();
    unsigned int getNumEigenFaces();
    unsigned int getNumTopEigenFaces();

    unsigned int getNumFacesAdded();
    unsigned int getNumFacesRemoved();

    Matrix getFaceEntry(int n); // Gets the nth flattened face in the face database
    Matrix getMeanFace();
    Matrix getEigenFaceEntry(int n); // Gets the nth flattened face in the eigen face database

    Matrix matchFace(Matrix& faceToMatch, double similarityScore); // faceToMatch must be a flattened face


    // todo display face, match face, remove face, and add face method

};
