// SONDA DE EXECUCAO EM OPENCL: a placa calcula o que o t8Flow precisaria?
//
// Nao basta o runtime estar instalado, e nao basta o `clinfo` listar o
// dispositivo: isso e' "presente", nao "funcional".  Esta sonda COMPILA um
// kernel, EXECUTA e CONFERE o resultado contra a CPU -- no mesmo padrao do
// `probe-c11.c` que decidiu a viabilidade do t8code.
//
// E o que ela calcula NAO e' um vetor-soma: e' o nucleo do LOTE DE ESTENCEIS que
// o usuario pediu.  Para cada ponto de consulta, com K pontos de suporte, monta a
// matriz normal do ajuste de minimos quadrados de ordem 1 em 2D (tres monomios:
// 1, dx, dy), resolve o 3x3 por eliminacao com pivotamento, e devolve o valor
// interpolado.  Mesma forma do problema real: muitos sistemas densos PEQUENOS.
//
// ORACULO, e e' o mesmo da suite: reproducao polinomial com CONTRASTE.  Um ajuste
// de ordem 1 reproduz exatamente campo LINEAR e erra no QUADRATICO.  Se os dois
// derem o mesmo erro, nao houve ajuste -- houve coincidencia, e a sonda diz isso.
//
// A PRECISAO E' O ACHADO QUE ESTA SONDA EXISTE PARA MEDIR.  Medido nesta maquina:
// a Intel UHD 630 tem `cl_khr_fp64`, e a Radeon via rusticl **NAO** -- ela nao
// expoe `CL_DEVICE_DOUBLE_FP_CONFIG`.  Entao a sonda roda em double onde ha' fp64
// e em float onde nao ha', e reporta o erro dos DOIS: a pergunta de projeto nao e'
// "a placa calcula", e' "a placa calcula com precisao suficiente para o MLS".
//
// Compilar:
//   cc -O2 bench/probe-opencl.c -o bench/probe-opencl -lOpenCL -lm
// Rodar (o rusticl precisa da variavel, senao a plataforma aparece com ZERO
// dispositivos -- que se le' como driver quebrado em vez de desabilitado):
//   RUSTICL_ENABLE=radeonsi ./bench/probe-opencl

#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>

// O <time.h> tem de vir com _POSIX_C_SOURCE para expor `clock_gettime` e
// `struct timespec`; definido aqui e nao na linha de compilacao para o arquivo
// compilar sozinho, com `cc probe-opencl.c -lOpenCL -lm`.
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif
#include <time.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define NCONS 20000   // pontos de consulta no lote
#define KSUP     16   // pontos de suporte por consulta
#define NPOLY     3   // 1, dx, dy  (ordem 1 em 2D)

// O kernel, em duas variantes de precisao.  `REAL` vem por -D na compilacao do
// programa OpenCL, nao por #ifdef no fonte: assim o MESMO texto serve aos dois
// casos e nao ha' chance de as variantes divergirem.
static const char *FONTE_KERNEL =
"#ifdef USA_DOUBLE\n"
"  #pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
"  typedef double real;\n"
"#else\n"
"  typedef float real;\n"
"#endif\n"
"\n"
"// Um estencil por item de trabalho: monta a matriz normal 3x3 do ajuste de\n"
"// ordem 1 sobre os KSUP pontos de suporte e resolve.  O valor interpolado em x\n"
"// e' o primeiro coeficiente, porque a base e' RELATIVA a x -- b(x) = (1,0,0).\n"
"// Base relativa nao e' detalhe de gosto: em coordenadas absolutas a matriz fica\n"
"// mal condicionada quando o dominio esta' longe da origem.\n"
"__kernel void estencil_lote(const int ksup,\n"
"                            __global const real *sup,   // 2*ksup por consulta\n"
"                            __global const real *val,   // ksup por consulta\n"
"                            __global const real *alvo,  // 2 por consulta\n"
"                            __global real *saida)\n"
"{\n"
"    const int q = get_global_id(0);\n"
"    const real ax = alvo[2*q], ay = alvo[2*q+1];\n"
"\n"
"    real M[3][3], r[3];\n"
"    for (int i = 0; i < 3; i++) { r[i] = 0; for (int j = 0; j < 3; j++) M[i][j] = 0; }\n"
"\n"
"    for (int k = 0; k < ksup; k++) {\n"
"        const real dx = sup[2*(q*ksup+k)] - ax;\n"
"        const real dy = sup[2*(q*ksup+k)+1] - ay;\n"
"        const real b[3] = { (real)1, dx, dy };\n"
"        const real f = val[q*ksup+k];\n"
"        for (int i = 0; i < 3; i++) {\n"
"            r[i] += b[i] * f;\n"
"            for (int j = 0; j < 3; j++) M[i][j] += b[i] * b[j];\n"
"        }\n"
"    }\n"
"\n"
"    // Gauss com pivotamento parcial.  Sem pivotamento, M[0][0] pode ser pequeno\n"
"    // quando o suporte e' simetrico em torno do alvo, e o resultado degrada em\n"
"    // silencio -- que e' o modo de falha que esta sonda existe para nao ter.\n"
"    for (int c = 0; c < 3; c++) {\n"
"        int piv = c;\n"
"        for (int i = c+1; i < 3; i++) if (fabs(M[i][c]) > fabs(M[piv][c])) piv = i;\n"
"        if (piv != c) {\n"
"            for (int j = 0; j < 3; j++) { real t = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = t; }\n"
"            real t = r[c]; r[c] = r[piv]; r[piv] = t;\n"
"        }\n"
"        for (int i = c+1; i < 3; i++) {\n"
"            const real m = M[i][c] / M[c][c];\n"
"            for (int j = c; j < 3; j++) M[i][j] -= m * M[c][j];\n"
"            r[i] -= m * r[c];\n"
"        }\n"
"    }\n"
"    real coef[3];\n"
"    for (int i = 2; i >= 0; i--) {\n"
"        real s = r[i];\n"
"        for (int j = i+1; j < 3; j++) s -= M[i][j] * coef[j];\n"
"        coef[i] = s / M[i][i];\n"
"    }\n"
"    saida[q] = coef[0];\n"
"}\n";

static double campo_linear(double x, double y) { return 1.0 + 2.0 * x - 3.0 * y; }
static double campo_quad(double x, double y)   { return campo_linear(x, y) + 4.0 * x * x; }

//! O MESMO calculo na CPU, em double -- o oraculo.
static double cpu_estencil(int ksup, const double *sup, const double *val, const double *alvo)
{
    double M[3][3] = {{0}}, r[3] = {0};
    for (int k = 0; k < ksup; k++) {
        const double dx = sup[2*k] - alvo[0], dy = sup[2*k+1] - alvo[1];
        const double b[3] = {1.0, dx, dy};
        for (int i = 0; i < 3; i++) {
            r[i] += b[i] * val[k];
            for (int j = 0; j < 3; j++) M[i][j] += b[i] * b[j];
        }
    }
    for (int c = 0; c < 3; c++) {
        int piv = c;
        for (int i = c+1; i < 3; i++) if (fabs(M[i][c]) > fabs(M[piv][c])) piv = i;
        if (piv != c) {
            for (int j = 0; j < 3; j++) { double t = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = t; }
            double t = r[c]; r[c] = r[piv]; r[piv] = t;
        }
        for (int i = c+1; i < 3; i++) {
            const double m = M[i][c] / M[c][c];
            for (int j = c; j < 3; j++) M[i][j] -= m * M[c][j];
            r[i] -= m * r[c];
        }
    }
    double coef[3];
    for (int i = 2; i >= 0; i--) {
        double s = r[i];
        for (int j = i+1; j < 3; j++) s -= M[i][j] * coef[j];
        coef[i] = s / M[i][i];
    }
    return coef[0];
}

static int tem_fp64(cl_device_id d)
{
    cl_device_fp_config fp = 0;
    clGetDeviceInfo(d, CL_DEVICE_DOUBLE_FP_CONFIG, sizeof fp, &fp, NULL);
    return fp != 0;
}

//! Roda o lote num dispositivo e devolve o maior erro contra a CPU, ou -1 se falhou.
static double roda(cl_device_id dev, int usa_double, int ksup,
                   const double *sup, const double *val, const double *alvo,
                   const double *esperado, int ncons, double *seg)
{
    cl_int e = CL_SUCCESS;
    cl_context ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &e);
    if (e != CL_SUCCESS) { printf("      clCreateContext falhou: %d\n", e); return -1; }
    cl_command_queue fila = clCreateCommandQueueWithProperties(ctx, dev, NULL, &e);
    if (e != CL_SUCCESS) { printf("      fila falhou: %d\n", e); clReleaseContext(ctx); return -1; }

    cl_program prog = clCreateProgramWithSource(ctx, 1, &FONTE_KERNEL, NULL, &e);
    const char *opcoes = usa_double ? "-DUSA_DOUBLE" : "";
    e = clBuildProgram(prog, 1, &dev, opcoes, NULL, NULL);
    if (e != CL_SUCCESS) {
        size_t n = 0;
        clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &n);
        char *log = malloc(n + 1);
        clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, n, log, NULL);
        log[n] = 0;
        printf("      COMPILACAO DO KERNEL FALHOU (%d):\n%s\n", e, log);
        free(log);
        clReleaseProgram(prog); clReleaseCommandQueue(fila); clReleaseContext(ctx);
        return -1;
    }
    cl_kernel k = clCreateKernel(prog, "estencil_lote", &e);

    const size_t nsup = (size_t)ncons * ksup * 2, nval = (size_t)ncons * ksup;
    const size_t tam = usa_double ? sizeof(double) : sizeof(float);
    void *hsup = malloc(nsup * tam), *hval = malloc(nval * tam), *halvo = malloc((size_t)ncons * 2 * tam);
    if (usa_double) {
        memcpy(hsup, sup, nsup * tam); memcpy(hval, val, nval * tam);
        memcpy(halvo, alvo, (size_t)ncons * 2 * tam);
    } else {
        for (size_t i = 0; i < nsup; i++) ((float*)hsup)[i] = (float)sup[i];
        for (size_t i = 0; i < nval; i++) ((float*)hval)[i] = (float)val[i];
        for (size_t i = 0; i < (size_t)ncons*2; i++) ((float*)halvo)[i] = (float)alvo[i];
    }

    cl_mem bsup = clCreateBuffer(ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, nsup*tam, hsup, &e);
    cl_mem bval = clCreateBuffer(ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, nval*tam, hval, &e);
    cl_mem balvo= clCreateBuffer(ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, (size_t)ncons*2*tam, halvo, &e);
    cl_mem bout = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY, (size_t)ncons*tam, NULL, &e);
    if (e != CL_SUCCESS) { printf("      buffers falharam: %d\n", e); return -1; }

    clSetKernelArg(k, 0, sizeof(int), &ksup);
    clSetKernelArg(k, 1, sizeof(cl_mem), &bsup);
    clSetKernelArg(k, 2, sizeof(cl_mem), &bval);
    clSetKernelArg(k, 3, sizeof(cl_mem), &balvo);
    clSetKernelArg(k, 4, sizeof(cl_mem), &bout);

    // Melhor de tres: a primeira execucao paga a compilacao em codigo de maquina
    // e o aquecimento da fila, e medir so' ela daria o custo de arranque como se
    // fosse o do calculo.
    double melhor = 1e300;
    for (int rep = 0; rep < 3; rep++) {
        size_t global = ncons;
        cl_event ev;
        clFinish(fila);
        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        e = clEnqueueNDRangeKernel(fila, k, 1, NULL, &global, NULL, 0, NULL, &ev);
        if (e != CL_SUCCESS) { printf("      enqueue falhou: %d\n", e); return -1; }
        clWaitForEvents(1, &ev);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        clReleaseEvent(ev);
        const double dt = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
        if (dt < melhor) melhor = dt;
    }
    *seg = melhor;

    void *hout = malloc((size_t)ncons * tam);
    e = clEnqueueReadBuffer(fila, bout, CL_TRUE, 0, (size_t)ncons*tam, hout, 0, NULL, NULL);
    if (e != CL_SUCCESS) { printf("      leitura falhou: %d\n", e); return -1; }

    double pior = 0.0;
    for (int q = 0; q < ncons; q++) {
        const double obtido = usa_double ? ((double*)hout)[q] : (double)((float*)hout)[q];
        const double d = fabs(obtido - esperado[q]);
        if (d > pior) pior = d;
    }

    free(hsup); free(hval); free(halvo); free(hout);
    clReleaseMemObject(bsup); clReleaseMemObject(bval);
    clReleaseMemObject(balvo); clReleaseMemObject(bout);
    clReleaseKernel(k); clReleaseProgram(prog);
    clReleaseCommandQueue(fila); clReleaseContext(ctx);
    return pior;
}

int main(void)
{
    // ---- o lote, gerado deterministicamente -------------------------------
    // Suporte em grade 4x4 ao redor do alvo, com o alvo FORA do centro de
    // simetria -- no centro, os termos impares se cancelam contra pesos
    // simetricos e um ajuste de ordem 1 reproduz quadratica por acidente.  Foi
    // medido na suite do t8Flow, e vale aqui igual.
    double *sup = malloc((size_t)NCONS * KSUP * 2 * sizeof(double));
    double *vlin = malloc((size_t)NCONS * KSUP * sizeof(double));
    double *vquad = malloc((size_t)NCONS * KSUP * sizeof(double));
    double *alvo = malloc((size_t)NCONS * 2 * sizeof(double));
    double *esp_lin = malloc((size_t)NCONS * sizeof(double));
    double *esp_quad = malloc((size_t)NCONS * sizeof(double));

    const double h = 0.1;
    for (int q = 0; q < NCONS; q++) {
        const double cx = 0.37 + 0.001 * (q % 211);
        const double cy = 0.44 + 0.001 * (q % 307);
        alvo[2*q] = cx; alvo[2*q+1] = cy;
        for (int k = 0; k < KSUP; k++) {
            const double px = cx + ((k % 4) - 1.3) * h;
            const double py = cy + ((k / 4) - 1.7) * h;
            sup[2*(q*KSUP+k)] = px;
            sup[2*(q*KSUP+k)+1] = py;
            vlin[q*KSUP+k] = campo_linear(px, py);
            vquad[q*KSUP+k] = campo_quad(px, py);
        }
        esp_lin[q] = campo_linear(cx, cy);
        esp_quad[q] = campo_quad(cx, cy);
    }

    // A CPU primeiro: ela e' o oraculo, e o seu proprio erro no campo linear diz
    // qual e' o piso de comparacao.
    double pior_cpu_lin = 0.0, pior_cpu_quad = 0.0;
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int q = 0; q < NCONS; q++) {
        const double v = cpu_estencil(KSUP, &sup[2*q*KSUP], &vlin[q*KSUP], &alvo[2*q]);
        const double d = fabs(v - esp_lin[q]);
        if (d > pior_cpu_lin) pior_cpu_lin = d;
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    const double seg_cpu = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
    for (int q = 0; q < NCONS; q++) {
        const double v = cpu_estencil(KSUP, &sup[2*q*KSUP], &vquad[q*KSUP], &alvo[2*q]);
        const double d = fabs(v - esp_quad[q]);
        if (d > pior_cpu_quad) pior_cpu_quad = d;
    }

    printf("lote: %d estenceis, %d pontos de suporte cada, ordem 1 em 2D\n", NCONS, KSUP);
    printf("\nCPU (double) -- o oraculo\n");
    printf("  campo linear     pior erro %.3e   %.4f s\n", pior_cpu_lin, seg_cpu);
    printf("  campo quadratico pior erro %.3e   <- CONTRASTE: ordem 1 nao reproduz grau 2\n",
           pior_cpu_quad);
    if (!(pior_cpu_quad > pior_cpu_lin * 1e3)) {
        printf("  AVISO: sem contraste na CPU (%.3e contra %.3e) -- o oraculo nao\n"
               "  discrimina, e nada abaixo significa nada\n", pior_cpu_lin, pior_cpu_quad);
    }

    // ---- os dispositivos ---------------------------------------------------
    cl_uint nplat = 0;
    clGetPlatformIDs(0, NULL, &nplat);
    if (nplat == 0) {
        printf("\nNENHUMA plataforma OpenCL.  O loader sozinho nao calcula: falta ICD em\n"
               "/etc/OpenCL/vendors/.  E o rusticl precisa de RUSTICL_ENABLE=radeonsi,\n"
               "senao aparece com ZERO dispositivos.\n");
        return 1;
    }
    cl_platform_id *plats = malloc(nplat * sizeof *plats);
    clGetPlatformIDs(nplat, plats, NULL);

    int achou = 0, falhou = 0;
    for (cl_uint p = 0; p < nplat; p++) {
        char pnome[256] = {0};
        clGetPlatformInfo(plats[p], CL_PLATFORM_NAME, sizeof pnome, pnome, NULL);
        cl_uint ndev = 0;
        clGetDeviceIDs(plats[p], CL_DEVICE_TYPE_ALL, 0, NULL, &ndev);
        if (ndev == 0) {
            printf("\nplataforma %s: ZERO dispositivos\n", pnome);
            continue;
        }
        cl_device_id *devs = malloc(ndev * sizeof *devs);
        clGetDeviceIDs(plats[p], CL_DEVICE_TYPE_ALL, ndev, devs, NULL);
        for (cl_uint d = 0; d < ndev; d++) {
            char dnome[256] = {0};
            clGetDeviceInfo(devs[d], CL_DEVICE_NAME, sizeof dnome, dnome, NULL);
            cl_uint cu = 0;
            clGetDeviceInfo(devs[d], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof cu, &cu, NULL);
            const int fp64 = tem_fp64(devs[d]);
            printf("\n%s | %s\n", pnome, dnome);
            printf("  %u unidades de calculo, fp64: %s\n", cu, fp64 ? "SIM" : "NAO");
            if (!fp64)
                printf("  sem fp64 o lote roda em FLOAT, e o erro abaixo diz se isso serve\n");

            achou++;
            double seg = 0.0;
            const double pior_lin = roda(devs[d], fp64, KSUP, sup, vlin, alvo, esp_lin, NCONS, &seg);
            if (pior_lin < 0) { falhou++; continue; }
            double seg2 = 0.0;
            const double pior_quad = roda(devs[d], fp64, KSUP, sup, vquad, alvo, esp_quad, NCONS, &seg2);

            printf("  campo linear     pior erro %.3e   %.4f s  (%.1fx a CPU)\n",
                   pior_lin, seg, seg > 0 ? seg_cpu / seg : 0.0);
            printf("  campo quadratico pior erro %.3e\n", pior_quad);
            if (pior_quad > pior_lin * 1e3)
                printf("  CONTRASTE ok: o ajuste esta' acontecendo na placa\n");
            else
                printf("  SEM CONTRASTE (%.3e contra %.3e): o resultado nao distingue\n"
                       "  ajuste de coincidencia -- nao confie neste dispositivo sem olhar\n",
                       pior_lin, pior_quad);
        }
        free(devs);
    }
    free(plats);

    // A RESSALVA QUE EVITA LER O GANHO COMO MAIOR DO QUE E'.  A CPU acima e'
    // SERIAL -- um nucleo, laco simples.  Numa maquina de N nucleos fisicos a
    // referencia honesta e' ~N vezes mais rapida, e o ganho da placa cai pelo
    // mesmo fator.  Dizer o ganho contra um nucleo sem dizer isso e' como medir
    // escalabilidade sem declarar quantos nucleos a maquina tem.
    {
        long nucleos = sysconf(_SC_NPROCESSORS_ONLN);
        printf("\nCOMO LER OS GANHOS: a CPU de referencia e' SERIAL (um nucleo).\n");
        printf("  esta maquina anuncia %ld cpus logicas -- e `nproc` conta THREADS, nao\n", nucleos);
        printf("  nucleos: confira com `lscpu`.  Com N nucleos a CPU paralela seria ~N vezes\n");
        printf("  mais rapida, e o ganho de cada placa cai pelo mesmo fator.  Um ganho de\n");
        printf("  14x contra um nucleo e' menos de 2x contra oito.\n");
        printf("\nE A PRECISAO E' DECISAO DE FORMULACAO, nao de engenharia: dispositivo sem\n");
        printf("  fp64 roda em float, e o erro do campo LINEAR acima e' o piso dele.  Se o\n");
        printf("  estencil alimenta o solver de pressao com tolerancia 1e-10, float nao\n");
        printf("  serve; se alimenta interpolacao de saida, serve.\n");
    }

    printf("\n%d dispositivo(s) exercitado(s), %d falharam\n", achou, falhou);
    return (achou > 0 && falhou == 0) ? 0 : 1;
}
