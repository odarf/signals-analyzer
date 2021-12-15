#include "inou.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>

namespace endian_io
{
  template <typename Word>
  std::ostream& write_word( std::ostream& outs, Word value, unsigned size = sizeof( Word ) )
  {
    for (; size; --size, value >>= 8)
      outs.put( static_cast <char> (value & 0xFF) );
    return outs;
  }
}

using namespace endian_io;

inou::inou(){}

QVector<double> inou::load(string path){
    QVector<double> outputData;
    ifstream file(path, ios::in | ios::binary);
    file.seekg(0, ios::end);
    int tam = file.tellg();
    tam /= 4;
    file.seekg(0, ios::beg);
    QVector<float> foo(tam);
    for (int i=0; i<tam; i++) {
        file.read(reinterpret_cast< char * >(&foo[i]), sizeof(float));
    }
    file.close();
    for (int i=0; i<foo.length(); i++) {
        outputData.append(static_cast<double>(foo[i]));
    }
    return outputData;
}

//Wav Header
struct wav_header_t
{
    char chunkID[4]; //"RIFF" = 0x46464952
    unsigned long chunkSize; //28 [+ sizeof(wExtraFormatBytes) + wExtraFormatBytes] + sum(sizeof(chunk.id) + sizeof(chunk.size) + chunk.size)
    char format[4]; //"WAVE" = 0x45564157
    char subchunk1ID[4]; //"fmt " = 0x20746D66
    unsigned long subchunk1Size; //16 [+ sizeof(wExtraFormatBytes) + wExtraFormatBytes]
    unsigned short audioFormat;
    unsigned short numChannels;
    unsigned long sampleRate;
    unsigned long byteRate;
    unsigned short blockAlign;
    unsigned short bitsPerSample;
    //[WORD wExtraFormatBytes;]
    //[Extra format bytes]
};

//Chunks
struct chunk_t
{
    char ID[4]; //"data" = 0x61746164
    unsigned long size;  //Chunk data bytes
};

typedef struct headder_file* header_p;

QVector<double> inou::loadWave(const char* fileName){
    QVector<double> outputData;
    FILE *fin = fopen(fileName, "rb");

    //Read WAV header
    wav_header_t header;
    fread(&header, sizeof(header), 1, fin);

    //Print WAV header
    std::clog << "WAV File Header read: " << std::endl;
    std::clog << "File Type: " << header.chunkID << std::endl;
    std::clog << "File Size: " <<  header.chunkSize << std::endl;
    std::clog << "WAV Marker: " <<  header.format << std::endl;
    std::clog << "Format Name: " <<  header.subchunk1ID << std::endl;
    std::clog << "Format Length: " <<  header.subchunk1Size << std::endl;
    std::clog << "Format Type: " <<  header.audioFormat << std::endl;
    std::clog << "Number of Channels: " <<  header.numChannels << std::endl;
    std::clog << "Sample Rate: " <<  header.sampleRate << std::endl;
    std::clog << "Sample Rate * Bits/Sample * Channels / 8: " <<  header.byteRate << std::endl;
    std::clog << "Bits per Sample * Channels / 8.1: " << header.blockAlign << std::endl;
    std::clog << "Bits per Sample: " << header.bitsPerSample << std::endl;

    //skip wExtraFormatBytes & extra format bytes
    //fseek(f, header.chunkSize - 16, SEEK_CUR);

    //Reading file
    chunk_t chunk;
    printf("id\t" "size\n");
    //go to data chunk
    while (true){
       fread(&chunk, sizeof(chunk), 1, fin);
       std::clog << chunk.ID[0] << chunk.ID[1] << chunk.ID[2] << chunk.ID[3] << " Chunk size: " << chunk.size << std::endl;
       if (*(unsigned int *)&chunk.ID == 0x61746164)
            break;
        //skip chunk data bytes
        fseek(fin, chunk.size, SEEK_CUR);
    }

    //Number of samples
    int sample_size = header.bitsPerSample / 8;
    int samples_count = (chunk.size * 8) / header.bitsPerSample;
    std::clog << "Samples count = " << samples_count << std::endl;

    short int *valueMain = new short int[samples_count];
    memset(valueMain, 0, sizeof(short int) * samples_count);

    //Reading data
    for (int i = 0; i < samples_count; i++){
        fread(&valueMain[i], sample_size, 1, fin);
    }

    fclose(fin);

    for (int i = 0; i < samples_count; i++){ //i < samples_count;
        outputData.append(valueMain[i]);
    }
    return outputData;
}

void inou::exportWave(QVector<double> inputdData, int samples_count, const char *fileToSave, double volume){
    std::ofstream fout(fileToSave, std::ios::binary);

    fout << "RIFF----WAVEfmt ";
    write_word(fout,     16, 4 );
    write_word(fout,      1, 2 );
    write_word(fout,      2, 2 );
    write_word(fout,  22050, 4 ); //22050,
    write_word(fout,  88200, 4 ); //88200
    write_word(fout,      4, 2 );
    write_word(fout,     16, 2 );

    //size_t data_chunk_pos = fout.tellp();
    fout << "data----";
    write_word(fout, samples_count);

    for (int i = 0; i < samples_count; i++){ //i < samples_count;
        write_word(fout, (int)(inputdData[i]*volume), 2); //Volume
    }

    /*size_t file_length = fout.tellp();
    fout.seekp(data_chunk_pos + 4);
    write_word(fout , file_length - data_chunk_pos + 8);
    fout.seekp(0 + 4);
    write_word(fout, file_length -8, 4);*/

    fout.close();
}
