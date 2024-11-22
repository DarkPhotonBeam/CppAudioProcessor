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

int main() {
    FILE *f;

    if (f = fopen("./sample.wav", "rb"); f == nullptr) {
        std::cerr << "Error opening file" << std::endl;
        return 1;
    }

    char chunkID[5];
    chunkID[4] = 0;
    unsigned int chunkSize;
    char format[5];
    format[4] = 0;

    fread(chunkID, 1, 4, f);
    fread(&chunkSize, 4, 1, f);
    fread(format, 1, 4, f);

    printf("Chunk ID: %s\n", chunkID);
    printf("File Size: %uB\n", chunkSize);
    printf("Format: %s\n", format);

    char subchunk1ID[5];
    subchunk1ID[4] = 0;
    unsigned int subchunk1Size;
    unsigned short audioFormat;
    unsigned short numChannels;
    unsigned int sampleRate;
    unsigned int byteRate;
    unsigned short blockAlign;
    unsigned short bitsPerSample;

    fread(subchunk1ID, 1, 4, f);
    fread(&subchunk1Size, 4, 1, f);
    fread(&audioFormat, 2, 1, f);
    fread(&numChannels, 2, 1, f);
    fread(&sampleRate, 4, 1, f);
    fread(&byteRate, 4, 1, f);
    fread(&blockAlign, 2, 1, f);
    fread(&bitsPerSample, 2, 1, f);

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
    char subchunk2ID[5];
    subchunk2ID[4] = 0;
    unsigned int subchunk2Size;

    fread(subchunk2ID, 1, 4, f);
    fread(&subchunk2Size, 4, 1, f);

    const size_t numSamples = subchunk2Size / blockAlign;

    printf("Subchunk2 ID: %s\n", subchunk2ID);
    printf("Subchunk2 Size: %uB\n", subchunk2Size);
    printf("Total number of samples: %lu\n", numSamples);



    unsigned short currL, currR;

    Eigen::Index offset = 10000;
    Eigen::Index fftBuffer = 1024 << 2;

    printf("Buffer Size: %ld\n", fftBuffer);

    Eigen::VectorXd lVec(fftBuffer);
    Eigen::VectorXd rVec(fftBuffer);

    for (Eigen::Index i = offset; i < offset + fftBuffer; ++i) {
        fread(&currL, bytesPerSample, 1, f);
        fread(&currR, bytesPerSample, 1, f);
        lVec(i-offset) = currL;
        rVec(i-offset) = currR;
    }

    printf("Read data.\n");

    Eigen::FFT<double> fft;
    printf("hmm\n");

    const Eigen::VectorXcd cl = fft.fwd(lVec);
    printf("Left FFT done.\n");
    const Eigen::VectorXcd cr = fft.fwd(rVec);
    printf("Right FFT done.\n");

    std::vector<double> vec{};

    int count = 0;
    int max_count = fftBuffer / 10;

    for (std::complex<double> z : cl) {
        if (count >= max_count) break;
        //std::cout << std::sqrt(z.imag()*z.imag() + z.real() * z.real()) << std::endl;
        vec.push_back(std::sqrt(z.imag()*z.imag() * z.real()*z.real()));
        ++count;
    }

    std::cout << "size: " << cl.size() << ", " << fftBuffer << std::endl;

    plt::plot(vec);
    plt::show();

    fclose(f);

    return 0;
}
