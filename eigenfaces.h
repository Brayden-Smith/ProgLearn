#include "data.h"

class FaceSpace {
private:
    Matrix faceDatabase;
    Matrix eigenFaces;
    long numFaces;
    unsigned int numEigenFaces;

    // For adding and removing faces (will come later)
    unsigned int numFacesAdded;
    unsigned int numFacesRemoved;

public:
    // Constructors
    FaceSpace(std::string& path);
    FaceSpace(Matrix faceMatrix);





};
