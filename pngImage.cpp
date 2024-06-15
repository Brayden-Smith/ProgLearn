#include "pngImage.h"

pngImage::pngImage(std::string filePath) {
    std::ifstream input(filePath, std::ios::binary);
    if(input) {
        input.seekg(0,std::ios::end);
        int size = input.tellg();
        input.seekg(0,std::ios::beg);

        //first 16 bytes are always the same, so we skip them
        input.seekg(16,std::ios::beg);

        //read in metadata :(
        char bytes[4];
        input.read(bytes,4);
        width = int((unsigned char)(bytes[0]) << 24 |
                    (unsigned char)(bytes[1]) << 16 |
                    (unsigned char)(bytes[2]) << 8 |
                    (unsigned char)(bytes[3]));
        input.read(bytes,4);
        height = int((unsigned char)(bytes[0]) << 24 |
                    (unsigned char)(bytes[1]) << 16 |
                    (unsigned char)(bytes[2]) << 8 |
                    (unsigned char)(bytes[3]));
        input.read(bytes,1);
        bitDepth = bytes[0];
        input.read(bytes,1);
        colorType = bytes[0];
        input.read(bytes,1);
        compressionMethod = bytes[0];
        input.read(bytes,1);
        filterMethod = bytes[0];
        input.read(bytes,1);
        interlaceMethod = bytes[0];

        int chunkLength = 0;
        std::string chunkType;

        bool iterate = true;
        while(iterate) {
            input.read(bytes,4);
            chunkLength = int((unsigned char)(bytes[0]) << 24 |
                        (unsigned char)(bytes[1]) << 16 |
                        (unsigned char)(bytes[2]) << 8 |
                        (unsigned char)(bytes[3]));

            input.read(bytes,4);
            chunkType = reinterpret_cast<char*>((bytes), sizeof(bytes));

            //if first letter is lower case than chunk is ancillary and not useful
            if(islower(chunkType[0])) {
                input.seekg(chunkLength + 4,std::ios::cur);
            }
            //end of file
            else if(chunkType == "IEND") {
                iterate = false;
            }
            else if(chunkType == "IDAT") {

            }
        }
        input.close();

    }

}