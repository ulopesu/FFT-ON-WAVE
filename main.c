
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#define GNUPLOT "gnuplot -persist"

typedef double complex cplx;

typedef struct {
    char ChunkID[4], Format[4], Subchunk1ID[4], Subchunk2ID[4];
    int ChunkSize, Subchunk1Size, SampleRate, ByteRate, Subchunk2Size;
    short AudioFormat, NumChannels, BlockAlign, BitsPerSample;
    int nSamples;
    short *Data;
} headerWAV;

short* createData(FILE *input, headerWAV header) {
    short *Data = (short*) malloc(header.nSamples * sizeof (short));
    return Data;
}

headerWAV readWAVHeader(FILE *input, headerWAV header) {
    //Leitura do arquivo WAV
    fread(header.ChunkID, 1, 4, input);
    fread(&header.ChunkSize, 4, 1, input);
    fread(header.Format, 1, 4, input);
    fread(header.Subchunk1ID, 1, 4, input);
    fread(&header.Subchunk1Size, 4, 1, input);
    fread(&header.AudioFormat, 2, 1, input);
    fread(&header.NumChannels, 2, 1, input);
    fread(&header.SampleRate, 4, 1, input);
    fread(&header.ByteRate, 4, 1, input);
    fread(&header.BlockAlign, 2, 1, input);
    fread(&header.BitsPerSample, 2, 1, input);
    fread(header.Subchunk2ID, 1, 4, input);
    fread(&header.Subchunk2Size, 4, 1, input);
    if ((header.ChunkSize - header.Subchunk2Size) > 36) {
        header.Subchunk2Size = header.ChunkSize - 36;
    } else {
        header.ChunkSize = header.Subchunk2Size + 36;
    }
    header.nSamples = (header.Subchunk2Size * 8) / (header.NumChannels * header.BitsPerSample);
    header.Data = createData(input, header);
    return header;
}

void readWavData(FILE *input, headerWAV header) {
    short sample;
    int i = 0;
    for (i = 0; i < header.nSamples; i++) {
        fread(&sample, header.BlockAlign, 1, input);
        header.Data[i] = sample;
    }
}

void normalize(cplx*v, int nSamples) {
    int i;
    int min = 0, max = 0;

    for (i = 0; i < nSamples; i++) {
        if (fabs(v[i]) < fabs(v[min])) min = i;
        if (fabs(v[i]) > fabs(v[max])) max = i;
    }

    double baseEscala = fabs(v[max]) - fabs(v[min]);
    double valor;

    for (i = 0; i < nSamples; i++) {
        valor = ((fabs(v[i]) - fabs(v[min])) / baseEscala);
        v[i] = fabs(valor);
    }
}

headerWAV writeWAV(FILE *output, headerWAV header) {
    fwrite(header.ChunkID, 1, 4, output);
    fwrite(&header.ChunkSize, 4, 1, output);
    fwrite(header.Format, 1, 4, output);
    fwrite(header.Subchunk1ID, 1, 4, output);
    fwrite(&header.Subchunk1Size, 4, 1, output);
    fwrite(&header.AudioFormat, 2, 1, output);
    fwrite(&header.NumChannels, 2, 1, output);
    fwrite(&header.SampleRate, 4, 1, output);
    fwrite(&header.ByteRate, 4, 1, output);
    fwrite(&header.BlockAlign, 2, 1, output);
    fwrite(&header.BitsPerSample, 2, 1, output);
    fwrite(header.Subchunk2ID, 1, 4, output);
    fwrite(&header.Subchunk2Size, 4, 1, output);
    fwrite(header.Data, header.BlockAlign, header.Subchunk2Size / header.BlockAlign, output);
    return header;
}

void fft(cplx* x, int N) {
    int estagios = log2(N), ka, indiceWD = 2, c, temp1, temp2;
    double w, s, t1, t2, temp;
    cplx c_temp, W;
    cplx * Xf = malloc(N * sizeof (cplx));
    int ant = 0, a, b, d;
    for (c = 1; c <= estagios; c++) {
        temp2 = N / c;
        for (temp1 = 0; temp1 < indiceWD; temp1++) {
            for (temp2 = 0; temp2 < (N / indiceWD); temp2++) {
                a = temp1 + (temp2 * indiceWD);
                d = indiceWD / 2;
                double angulo = temp1 * 2 * M_PI / indiceWD;
                W = cos(angulo) -(sin(angulo) * I);
                if (temp1 >= d) {
                    b = a - d;
                    Xf[a] = x[b] + x[a] * W;
                } else {
                    b = a + d;
                    Xf[a] = x[a] + x[b] * W;
                }
            }
        }
        //ATUALIZA AS ENTRADAS P/ O PROXIMO ESTAGIO
        for (ka = 0; ka < N; ka++) {
            x[ka] = Xf[ka];
        }
        indiceWD = indiceWD * 2;
    }
}

void ifft(cplx* x, int N) {
    int estagios = log2(N), ka, indiceWD = 2, c, temp1, temp2;
    double w, s, t1, t2, temp, angulo;
    cplx c_temp, W;
    cplx * Xf = malloc(N * sizeof (cplx));
    int ant = 0, a, b, d;
    for (c = 1; c <= estagios; c++) {
        temp2 = N / c;
        angulo = temp1 * 2 * M_PI / indiceWD;
        W = cos(angulo) + (sin(angulo) * I);
        for (temp1 = 0; temp1 < indiceWD; temp1++) {
            for (temp2 = 0; temp2 < (N / indiceWD); temp2++) {
                a = temp1 + (temp2 * indiceWD);
                d = indiceWD / 2;
                if (temp1 >= d) {
                    b = a - d;
                    Xf[a] = x[b] + x[a] * W;
                } else {
                    b = a + d;
                    Xf[a] = x[a] + x[b] * W;
                }
            }
        }
        //ATUALIZA AS ENTRADAS P/ O PROXIMO ESTAGIO
        for (ka = 0; ka < N; ka++) {
            x[ka] = Xf[ka];
        }
        indiceWD = indiceWD * 2;
    }
}

void show(const char * s, cplx *buf, int n) {
    printf("%s", s);
    for (int i = 0; i < n; i++)
        if (!cimag(buf[i]))
            printf("%g ", creal(buf[i]));
        else
            printf("(%g, %g) ", creal(buf[i]), cimag(buf[i]));
}

int main(int argc, char **argv) {
    int i, length;
    headerWAV header, headerOut;
    cplx *v;

    FILE *input = fopen(argv[1], "rb"), *output = fopen("saida.wav", "wb");

    if (input == NULL)
        printf("Arquivo nao pode ser aberto!!!\n");

    header = readWAVHeader(input, header);

    if (output == NULL)
        printf("Arquivo nao pode ser aberto!!!\n");

    readWavData(input, header);

    v = (cplx*) malloc(sizeof (cplx) * header.nSamples);
    for (i = 0; i < header.nSamples; i++) {
        v[i] = header.Data[i];
    }


    FILE *dados0 = fopen("dados0.txt", "wb"), *dados1 = fopen("dados1.txt", "wb"), *dados2 = fopen("dados2.txt", "wb");

    fprintf(dados0, "# x \t f(x) \n");
    for (i = 0; i < header.nSamples; i++) {
        fprintf(dados0, "%d \t %f \n", i, creal(v[i]));
    }

    //normalize(v, header.nSamples);
    fft(v, header.nSamples);

    fprintf(dados1, "# x \t f(x) \n");
    for (i = 0; i < header.nSamples; i++) {
        fprintf(dados1, "%f \t %f \n", creal(v[i]), cimag(v[i]));
    }

    fprintf(dados2, "# x \t f(x) \n");
    double x = (header.SampleRate / (double) header.nSamples), y;
    for (i = 0; i < header.nSamples; i++) {
        y = x * (i + 1);
        fprintf(dados2, "%f \t %f \n", y, fabs(v[i]));
    }

    /*
    ifft(v, header.nSamples);
    for (i = 0; i < header.nSamples; i++) {
        header.Data[i] = v[i];
    }
     */

    headerOut = writeWAV(output, header);

    fclose(dados0);
    fclose(dados1);
    fclose(dados2);
    fclose(input);
    fclose(output);

    FILE *gp0, *gp1, *gp2;
    //char setXRange[100];

    gp0 = popen(GNUPLOT, "w");
    gp1 = popen(GNUPLOT, "w");
    gp2 = popen(GNUPLOT, "w");

    if (gp0 == NULL || gp1 == NULL || gp2 == NULL) {
        printf("Erro ao abrir pipe para o GNU plot.\n"
                "Instale com 'sudo apt-get install gnuplot'\n");
        exit(0);
    }
    /*
    char * commands0[] = {"set title \"HISTOGRAMA DOS DADOS NO TEMPO (ANTES DA FFT)\"", "set style data histogram", "plot 'dados0.txt' using 2"};
    for (i = 0; i < 3; i++) {
        fprintf(gp0, "%s \n", commands0[i]);
    }

    char * commands1[] = {"set title \"PARTE REAL E IMAGINARIA (APOS FFT)\"", "plot 'dados1.txt'"};
    for (i = 0; i < 2; i++) {
        fprintf(gp1, "%s \n", commands1[i]);
    }
     */
    //sprintf(setXRange, "set xrange [0:%d+70000]", (int) (header.SampleRate));

    char * commands2[] = {"set title \"AMPLITUDE X FREQUENCIA - (APOS FFT)\"", "plot 'dados2.txt'"};
    for (i = 0; i < 2; i++) {
        fprintf(gp2, "%s \n", commands2[i]);
    }


    fclose(gp0);
    fclose(gp1);
    fclose(gp2);
    return 0;
}