#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

unsigned int switchEndianness(unsigned int x) {
    unsigned int byte0 = (x & 0x000000ff) << 24;
    unsigned int byte1 = (x & 0x0000ff00) << 8;
    unsigned int byte2 = (x & 0x00ff0000) >> 8;
    unsigned int byte3 = (x & 0xff000000) >> 24;
    return byte3 | byte2 | byte1 | byte0;
}

size_t freadwrite(
    void *__ptr,
    size_t __size,
    size_t __nitems,
    FILE *__srcStream, FILE *__targetStream) {

    size_t bytesRead = fread(__ptr, __size, __nitems, __srcStream);
    size_t bytesWritten = fwrite(__ptr, __size, __nitems, __targetStream);

    return bytesRead + bytesWritten;
}

int main() {
    FILE *f;
    FILE *dest;

    if (f = fopen("./audio/WhiteNoise.wav", "rb"); f == nullptr) {
        std::cerr << "Error opening file" << std::endl;
        return 1;
    }

    if (dest = fopen("./filtered.wav", "wb"); f == nullptr) {
        std::cerr << "Error opening dest file" << std::endl;
        return 1;
    }

    char chunkID[5];
    chunkID[4] = 0;
    unsigned int chunkSize;
    char format[5];
    format[4] = 0;

    int bruh[17];

    //

    freadwrite(chunkID, 1, 4, f, dest);
    freadwrite(&chunkSize, 4, 1, f, dest);
    freadwrite(format, 1, 4, f, dest);

    printf("Chunk ID: %s\n", chunkID);
    printf("File Size: %uB\n", chunkSize);
    printf("Format: %s\n", format);

    char subchunk0ID[5];
    subchunk0ID[4] = 0;
    unsigned int subchunk0Size;

    freadwrite(subchunk0ID, 1, 4, f, dest);
    freadwrite(&subchunk0Size, 4, 1, f, dest);

    printf("Subchunk0 ID: %s\n", subchunk0ID);
    printf("Subchunk0 Size: %uB\n", subchunk0Size);

    // TODO: handle both types of .wav files

    freadwrite(bruh, 4, 16, f, dest);


    //fread(bruh, 4, 5, f);

    char subchunk1ID[5];
    subchunk1ID[4] = 0;
    unsigned int subchunk1Size;
    unsigned short audioFormat;
    unsigned short numChannels;
    unsigned int sampleRate;
    unsigned int byteRate;
    unsigned short blockAlign;
    unsigned short bitsPerSample;

    freadwrite(subchunk1ID, 1, 4, f, dest);
    freadwrite(&subchunk1Size, 4, 1, f, dest);
    freadwrite(&audioFormat, 2, 1, f, dest);
    freadwrite(&numChannels, 2, 1, f, dest);
    freadwrite(&sampleRate, 4, 1, f, dest);
    freadwrite(&byteRate, 4, 1, f, dest);
    freadwrite(&blockAlign, 2, 1, f, dest);
    freadwrite(&bitsPerSample, 2, 1, f, dest);

    printf("Subchunk1 ID: %s\n", subchunk1ID);
    printf("Subchunk1 Size: %uB\n", subchunk1Size);
    printf("Audio format: %hu (%s)\n", audioFormat, (audioFormat == 1 ? "PCM" : "Non-PCM"));
    printf("Num Channels: %hu\n", numChannels);
    printf("Sample Rate: %u\n", sampleRate);
    printf("Byte Rate: %u\n", byteRate);
    printf("Bit Rate: %u kb/s\n", byteRate * 8 / 1000);
    printf("Block Align: %hu\n", blockAlign);
    printf("Bits Per Sample: %hu\n", bitsPerSample);

    unsigned short bytesPerSample = bitsPerSample / 8;
    printf("Bytes per Sample: %hu\n", bytesPerSample);


    char subchunk2ID[5];
    subchunk2ID[4] = 0;
    unsigned int subchunk2Size;

    freadwrite(subchunk2ID, 1, 4, f, dest);
    freadwrite(&subchunk2Size, 4, 1, f, dest);

    const size_t numSamples = subchunk2Size / blockAlign;

    printf("Subchunk2 ID: %s\n", subchunk2ID);
    printf("Subchunk2 Size: %uB\n", subchunk2Size);
    printf("Total number of samples: %lu\n", numSamples);

    std::vector<unsigned int> dataLeft;
    dataLeft.reserve(numSamples);
    std::vector<unsigned int> dataRight;
    dataRight.reserve(numSamples);


    for (unsigned int i = 0; i < numSamples; ++i) {
        //unsigned short currL, currR;
        fread(&dataLeft[i], bytesPerSample, 1, f);
        fread(&dataRight[i], bytesPerSample, 1, f);

        if (i % 100000 == 0) {
            printf("%u read %u\n", i, dataLeft[i]);
        }

        //fwrite(&currL, bytesPerSample, 1, dest);
        //fwrite(&currR, bytesPerSample, 1, dest);
    }

    size_t samplesLeft = numSamples;

    Eigen::Index offset = 10000;
    //long bufsz = 1 << 10;
    Eigen::Index fftBuffer = 1 << 15;
    printf("Buffer Size: %ld\n", fftBuffer);

    Eigen::FFT<double> fft;



    while (samplesLeft > 0) {
        if (samplesLeft < fftBuffer) {
            fftBuffer = samplesLeft;
        }
        printf("%lu samples left\n", samplesLeft);

        const double buf = static_cast<double>(fftBuffer);
        const double ctf = buf / 2.1;

        const int cutoff = static_cast<int>(ctf);
        const int plotCutoff = static_cast<int>(ctf);

        size_t begin = numSamples - samplesLeft;
        samplesLeft -= fftBuffer;

        Eigen::VectorXd lVec(fftBuffer);
        Eigen::VectorXd rVec(fftBuffer);

        for (Eigen::Index i = 0; i < fftBuffer; ++i) {
            lVec(i) = dataLeft[begin+i];
            rVec(i) = dataRight[begin+i];
        }

        Eigen::VectorXcd cl = fft.fwd(lVec);
        Eigen::VectorXcd cr = fft.fwd(rVec);

        std::vector<double> plotData{};

        const Eigen::VectorXd::Index m = lVec.size() / 2;

        for (Eigen::Index i = 0; i < cl.size(); ++i) {
            if (i < plotCutoff) plotData.push_back(std::sqrt(cl(i).real()*cl(i).real()+cr(i).imag()*cr(i).imag()));
        }

        for (int j = -cutoff; j <= cutoff; ++j) {
            if (m + j < 0 || m + j >= fftBuffer) continue;
            cl(m + j) = std::complex<double>(0, 0);
            cr(m + j) = std::complex<double>(0, 0);
        }

        Eigen::VectorXd invL = fft.inv(cl).real();
        Eigen::VectorXd invR = fft.inv(cr).real();

        for (Eigen::Index i = 0; i < fftBuffer; ++i) {
            dataLeft[begin+i] = static_cast<unsigned int>(invL(i));
            dataRight[begin+i] = static_cast<unsigned int>(invR(i));
        }

        std::string s =  "./plots/plot";
        s.append(std::to_string(samplesLeft));
        s.append(".png");

        plt::plot(plotData);
        plt::save(s);
    }

    // write back
    for (unsigned int i = 0; i < numSamples; ++i) {
        //unsigned short l = dataLeft[i];
        //unsigned short r = dataRight[i];
        if (i % 100000 == 0) {
            printf("%u write %u\n", i, dataLeft[i]);
        }
        fwrite(&dataLeft[i], bytesPerSample, 1, dest);
        fwrite(&dataRight[i], bytesPerSample, 1, dest);
    }

    fclose(f);
    fclose(dest);

    return 0;
}
