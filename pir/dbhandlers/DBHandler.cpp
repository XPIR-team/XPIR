#include "DBHandler.hpp"

void DBHandler::readAggregatedStream(uint64_t streamNb, uint64_t alpha, uint64_t offset, uint64_t bytes_per_file, char* rawBits) {
    uint64_t fileByteSize = std::min(bytes_per_file, getmaxFileBytesize()-offset);
    uint64_t startStream = streamNb*alpha;
    uint64_t endStream = std::min(streamNb*alpha + alpha - 1, getNbStream() - 1);
    uint64_t paddingStreams = (streamNb*alpha+alpha) >= getNbStream() ? (streamNb*alpha+alpha) - getNbStream() : 0;

  #pragma omp critical
    {
        for (int i=startStream; i <= endStream; i++)
        {
            openStream(i, offset);

            // Just read the file (plus padding for that file)
            readStream(i, rawBits + (i % alpha) * fileByteSize, fileByteSize);

            closeStream(i);
        }

        if(paddingStreams !=0)
        {
            bzero(rawBits + (endStream % alpha) * fileByteSize, fileByteSize*paddingStreams);
        }
    }
}
