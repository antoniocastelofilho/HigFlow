// Implementacao do lote.  O PORQUE da divisao CPU/GPU e da forma do resultado
// esta' no `lote2d.hxx`.  Aqui ficam as decisoes de implementacao:
//
//   1. O BACK-END DE CPU E' O NUCLEO HERDADO em laco -- `Mls::pesos` por consulta.
//      Nao ha' uma segunda implementacao do ajuste na CPU, e isso e' deliberado:
//      duas implementacoes do mesmo calculo divergem, e aqui a CPU e' o oraculo.
//      Se ela fosse um codigo novo, nao haveria oraculo nenhum.
//   2. O KERNEL DE GPU resolve o sistema normal em double quando o dispositivo
//      tem fp64, e em float quando nao tem.  Medido nesta maquina: a Intel UHD 630
//      tem fp64 e reproduz identico a' CPU (6,66e-16); a Radeon via rusticl NAO
//      tem, roda em float, e o erro sobe para 2,8e-07.  Por isso a precisao do
//      dispositivo e' REPORTADA e nao escondida.
//   3. TODO O CODIGO DE OPENCL ESTA' ATRAS DE `T8FLOW_COM_OPENCL`.  Sem a macro,
//      este arquivo compila sem incluir `CL/cl.h` e a classe funciona so' com CPU
//      -- que e' o caso de quem nao tem placa, e e' o caso que tem de continuar
//      funcionando.
#include "lote2d.hxx"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#ifdef T8FLOW_COM_OPENCL
#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>
// O EPSMACH do nucleo herdado -- o kernel TEM de usar o mesmo, senao os dois
// back-ends calculam ajustes diferentes.
#include "utils.h"
#endif

namespace t8flow2d {

int LoteEstenceis::num_monomios(int ordem)
{
  // Em 2D: (ordem+1)(ordem+2)/2 -- 3 para ordem 1, 6 para ordem 2, 10 para 3.
  return (ordem + 1) * (ordem + 2) / 2;
}

// ---------------------------------------------------------------------------
// O lado da GPU, inteiro atras da macro
// ---------------------------------------------------------------------------

#ifdef T8FLOW_COM_OPENCL

namespace {

//! O kernel: um estencil por item de trabalho.
//
// A BASE E' RELATIVA AO ALVO, e isso nao e' gosto: em coordenadas absolutas a
// matriz normal fica mal condicionada quando o dominio esta' longe da origem.
// Com a base relativa, o valor interpolado em x e' o primeiro coeficiente, porque
// b(x) = (1, 0, 0, ...).
//
// NMON entra por -D, nao por #if no texto: o mesmo fonte serve as duas ordens, e
// assim as variantes nao podem divergir.
const char *FONTE =
"#ifdef USA_DOUBLE\n"
"  #pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
"  typedef double real;\n"
"#else\n"
"  typedef float real;\n"
"#endif\n"
"\n"
"static void base(real dx, real dy, real *b)\n"
"{\n"
"    b[0] = (real)1; b[1] = dx; b[2] = dy;\n"
"#if NMON > 3\n"
"    b[3] = dx*dx; b[4] = dx*dy; b[5] = dy*dy;\n"
"#endif\n"
"}\n"
"\n"
"// Pesos do ajuste PONDERADO, e a ponderacao NAO e' opcional: o `wls` herdado e'\n"
"// WLSMOVING, com peso D_i = 1/sqrt(|x-p_i|^2 + EPSMACH) aplicado as linhas.  Com\n"
"// 16 pontos e 3 monomios o sistema tem 13 graus de liberdade, entao mudar o peso\n"
"// muda os PESOS inteiros -- e os dois ajustes continuam reproduzindo os mesmos\n"
"// polinomios, o que torna a reproducao incapaz de acusar a divergencia.  Quem\n"
"// acusou foi a comparacao GPU-contra-CPU: 0,311 de diferenca.\n"
"//\n"
"// w_k = D_k^2 * e0^T M^{-1} b_k,  com M = sum_i D_i^2 b_i b_i^T.  Resolve-se UM\n"
"// sistema por consulta e aplica-se a cada ponto, em vez de montar a inversa.\n"
"__kernel void lote(const int ksup,\n"
"                   __global const real *sup,\n"
"                   __global const real *alvo,\n"
"                   __global real *pesos)\n"
"{\n"
"    const int q = get_global_id(0);\n"
"    const real ax = alvo[2*q], ay = alvo[2*q+1];\n"
"\n"
"    real M[NMON][NMON], b[NMON];\n"
"    for (int i = 0; i < NMON; i++) for (int j = 0; j < NMON; j++) M[i][j] = 0;\n"
"    for (int k = 0; k < ksup; k++) {\n"
"        const real dx = sup[2*(q*ksup+k)] - ax;\n"
"        const real dy = sup[2*(q*ksup+k)+1] - ay;\n"
"        base(dx, dy, b);\n"
"        const real p2 = (real)1 / (dx*dx + dy*dy + (real)EPSMACH);   // = D_k^2\n"
"        for (int i = 0; i < NMON; i++)\n"
"            for (int j = 0; j < NMON; j++) M[i][j] += p2 * b[i] * b[j];\n"
"    }\n"
"\n"
"    // Resolve M y = e0 por Gauss com pivotamento parcial.  Sem pivotamento,\n"
"    // M[0][0] pode ser pequeno com suporte simetrico e o resultado degrada em\n"
"    // silencio.\n"
"    real y[NMON];\n"
"    for (int i = 0; i < NMON; i++) y[i] = (i == 0) ? (real)1 : (real)0;\n"
"    for (int c = 0; c < NMON; c++) {\n"
"        int piv = c;\n"
"        for (int i = c+1; i < NMON; i++) if (fabs(M[i][c]) > fabs(M[piv][c])) piv = i;\n"
"        if (piv != c) {\n"
"            for (int j = 0; j < NMON; j++) { real t = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = t; }\n"
"            real t = y[c]; y[c] = y[piv]; y[piv] = t;\n"
"        }\n"
"        for (int i = c+1; i < NMON; i++) {\n"
"            const real m = M[i][c] / M[c][c];\n"
"            for (int j = c; j < NMON; j++) M[i][j] -= m * M[c][j];\n"
"            y[i] -= m * y[c];\n"
"        }\n"
"    }\n"
"    for (int i = NMON-1; i >= 0; i--) {\n"
"        real s = y[i];\n"
"        for (int j = i+1; j < NMON; j++) s -= M[i][j] * y[j];\n"
"        y[i] = s / M[i][i];\n"
"    }\n"
"\n"
"    for (int k = 0; k < ksup; k++) {\n"
"        const real dx = sup[2*(q*ksup+k)] - ax;\n"
"        const real dy = sup[2*(q*ksup+k)+1] - ay;\n"
"        base(dx, dy, b);\n"
"        const real p2 = (real)1 / (dx*dx + dy*dy + (real)EPSMACH);\n"
"        real w = 0;\n"
"        for (int i = 0; i < NMON; i++) w += y[i] * b[i];\n"
"        pesos[q*ksup+k] = p2 * w;\n"
"    }\n"
"}\n";

}  // namespace

//! O estado de OpenCL, criado uma vez.  Falha em qualquer etapa deixa o ponteiro
//! nulo e a classe inteira cai na CPU -- sem excecao, porque nao ter GPU e' um
//! estado NORMAL e nao um erro.
struct LoteEstenceis::Gpu {
  cl_context ctx = nullptr;
  cl_command_queue fila = nullptr;
  cl_program prog = nullptr;
  cl_kernel kernel = nullptr;
  cl_device_id dev = nullptr;
  bool fp64 = false;
  std::string nome;
  std::string porque_nao;
};

#else   // sem OpenCL em tempo de compilacao

struct LoteEstenceis::Gpu {
  std::string porque_nao = "compilado SEM OpenCL (defina T8FLOW_COM_OPENCL)";
};

#endif

// ---------------------------------------------------------------------------
// LoteEstenceis
// ---------------------------------------------------------------------------

LoteEstenceis::LoteEstenceis(int ordem, int max_pontos)
    : ordem_(ordem), max_pontos_(max_pontos), mls_(ordem, max_pontos)
{
#ifdef T8FLOW_COM_OPENCL
  gpu_ = new Gpu();

  // A GPU so' atende as ordens que o kernel tem base para: 1 e 2 em 2D.  Ordem
  // maior NAO cai na GPU em silencio -- `gpu_disponivel()` devolve falso e o
  // diagnostico diz o motivo.
  const int nmon = num_monomios(ordem);
  if (nmon != 3 && nmon != 6) {
    gpu_->porque_nao = "o kernel cobre ordem 1 e 2 em 2D (3 ou 6 monomios); esta e' " +
                       std::to_string(nmon);
    return;
  }

  cl_uint nplat = 0;
  if (clGetPlatformIDs(0, nullptr, &nplat) != CL_SUCCESS || nplat == 0) {
    gpu_->porque_nao = "nenhuma plataforma OpenCL -- falta ICD em /etc/OpenCL/vendors/"
                       " (e o rusticl precisa de RUSTICL_ENABLE=radeonsi)";
    return;
  }
  std::vector<cl_platform_id> plats(nplat);
  clGetPlatformIDs(nplat, plats.data(), nullptr);

  // Prefere dispositivo com fp64: a precisao vale mais que a velocidade quando o
  // resultado alimenta o solver.  Medido: a placa mais rapida desta maquina e'
  // justamente a que NAO tem fp64.
  cl_device_id escolhido = nullptr;
  bool escolhido_fp64 = false;
  std::string escolhido_nome;
  for (cl_uint p = 0; p < nplat; p++) {
    cl_uint ndev = 0;
    if (clGetDeviceIDs(plats[p], CL_DEVICE_TYPE_ALL, 0, nullptr, &ndev) != CL_SUCCESS || ndev == 0)
      continue;
    std::vector<cl_device_id> devs(ndev);
    clGetDeviceIDs(plats[p], CL_DEVICE_TYPE_ALL, ndev, devs.data(), nullptr);
    for (cl_uint d = 0; d < ndev; d++) {
      cl_device_fp_config fp = 0;
      clGetDeviceInfo(devs[d], CL_DEVICE_DOUBLE_FP_CONFIG, sizeof fp, &fp, nullptr);
      char nome[256] = {0};
      clGetDeviceInfo(devs[d], CL_DEVICE_NAME, sizeof nome, nome, nullptr);
      const bool tem = (fp != 0);
      if (escolhido == nullptr || (tem && !escolhido_fp64)) {
        escolhido = devs[d];
        escolhido_fp64 = tem;
        escolhido_nome = nome;
      }
    }
  }
  if (escolhido == nullptr) {
    gpu_->porque_nao = "plataforma OpenCL presente, porem ZERO dispositivos";
    return;
  }

  cl_int e = CL_SUCCESS;
  gpu_->ctx = clCreateContext(nullptr, 1, &escolhido, nullptr, nullptr, &e);
  if (e != CL_SUCCESS) { gpu_->porque_nao = "clCreateContext falhou"; return; }
  gpu_->fila = clCreateCommandQueueWithProperties(gpu_->ctx, escolhido, nullptr, &e);
  if (e != CL_SUCCESS) { gpu_->porque_nao = "fila de comandos falhou"; return; }

  // EPSMACH vem do herdado (`utils.h:46`), por -D e nao copiado no texto do
  // kernel: um numero magico repetido em dois lugares e' um numero magico que vai
  // divergir.
  const std::string opcoes = std::string(escolhido_fp64 ? "-DUSA_DOUBLE " : "") +
                             "-DNMON=" + std::to_string(nmon) +
                             " -DEPSMACH=" + std::to_string(EPSMACH);
  gpu_->prog = clCreateProgramWithSource(gpu_->ctx, 1, &FONTE, nullptr, &e);
  e = clBuildProgram(gpu_->prog, 1, &escolhido, opcoes.c_str(), nullptr, nullptr);
  if (e != CL_SUCCESS) {
    size_t n = 0;
    clGetProgramBuildInfo(gpu_->prog, escolhido, CL_PROGRAM_BUILD_LOG, 0, nullptr, &n);
    std::string log(n, '\0');
    clGetProgramBuildInfo(gpu_->prog, escolhido, CL_PROGRAM_BUILD_LOG, n, &log[0], nullptr);
    gpu_->porque_nao = "o kernel nao compilou: " + log;
    return;
  }
  gpu_->kernel = clCreateKernel(gpu_->prog, "lote", &e);
  if (e != CL_SUCCESS) { gpu_->porque_nao = "clCreateKernel falhou"; return; }

  gpu_->dev = escolhido;
  gpu_->fp64 = escolhido_fp64;
  gpu_->nome = escolhido_nome + (escolhido_fp64 ? " (fp64)" : " (float -- SEM fp64)");

  // QUEM COMPILOU ESCOLHEU: com `make GPU=sim` o lote JA' comeca em GPU, sem o
  // chamador pedir.  Nos outros modos a placa fica disponivel e o padrao e' CPU.
#ifdef T8FLOW_GPU_PADRAO
  back_end_ = BackEnd::Gpu;
#endif
#else
  gpu_ = new Gpu();
#endif

  // PEDIRAM GPU E NAO HA': avisa em stderr e segue em CPU.  Nao lanca, porque
  // derrubar o programa por falta de placa seria pior que calcular mais devagar --
  // mas nao cala, porque a escolha de quem compilou foi ignorada e isso tem de
  // aparecer.  O aviso e' uma vez por objeto, no arranque, e nao no laco quente.
  if (gpu_pedida_na_compilacao() && !gpu_disponivel()) {
    std::fprintf(stderr,
                 "t8flow: compilado com GPU=sim, mas o lote vai rodar em CPU -- %s\n",
                 diagnostico_gpu().c_str());
  }
}

LoteEstenceis::~LoteEstenceis()
{
#ifdef T8FLOW_COM_OPENCL
  if (gpu_ != nullptr) {
    if (gpu_->kernel) clReleaseKernel(gpu_->kernel);
    if (gpu_->prog) clReleaseProgram(gpu_->prog);
    if (gpu_->fila) clReleaseCommandQueue(gpu_->fila);
    if (gpu_->ctx) clReleaseContext(gpu_->ctx);
  }
#endif
  delete gpu_;
}

bool LoteEstenceis::gpu_pedida_na_compilacao()
{
#ifdef T8FLOW_GPU_PADRAO
  return true;
#else
  return false;
#endif
}

bool LoteEstenceis::gpu_disponivel() const
{
#ifdef T8FLOW_COM_OPENCL
  return gpu_ != nullptr && gpu_->kernel != nullptr;
#else
  return false;
#endif
}

std::string LoteEstenceis::diagnostico_gpu() const
{
  if (gpu_disponivel()) {
#ifdef T8FLOW_COM_OPENCL
    return gpu_->nome;
#endif
  }
  return gpu_ != nullptr ? gpu_->porque_nao : "estado de GPU ausente";
}

void LoteEstenceis::usa_back_end(BackEnd b)
{
  if (b == BackEnd::Gpu && !gpu_disponivel())
    throw std::runtime_error("LoteEstenceis: GPU pedida mas indisponivel -- " +
                             diagnostico_gpu());
  back_end_ = b;
}

void LoteEstenceis::calcula(const std::vector<ConsultaEstencil> &consultas,
                            std::vector<Estencil> &saida) const
{
  // A VALIDACAO E' A MESMA NOS DOIS BACK-ENDS, e vem antes de escolher: suporte
  // invalido tem de ser recusado igual, senao "funciona na CPU e quebra na GPU"
  // passa a ser possivel.
  const int minimo = Mls::min_pontos(ordem_);
  for (size_t q = 0; q < consultas.size(); q++) {
    const ConsultaEstencil &c = consultas[q];
    if (c.ids.size() != c.suporte.size())
      throw std::invalid_argument("LoteEstenceis: consulta " + std::to_string(q) + " tem " +
                                  std::to_string(c.ids.size()) + " ids para " +
                                  std::to_string(c.suporte.size()) + " pontos");
    if (static_cast<int>(c.suporte.size()) < minimo)
      throw std::invalid_argument("LoteEstenceis: consulta " + std::to_string(q) +
                                  " com suporte de " + std::to_string(c.suporte.size()) +
                                  " pontos, minimo " + std::to_string(minimo) +
                                  " para ordem " + std::to_string(ordem_));
    if (static_cast<int>(c.suporte.size()) > max_pontos_)
      throw std::invalid_argument("LoteEstenceis: consulta " + std::to_string(q) +
                                  " excede max_pontos " + std::to_string(max_pontos_));
  }

  saida.assign(consultas.size(), Estencil{});
  if (consultas.empty()) return;

  if (back_end_ == BackEnd::Gpu && gpu_disponivel()) calcula_gpu(consultas, saida);
  else calcula_cpu(consultas, saida);
}

void LoteEstenceis::calcula_cpu(const std::vector<ConsultaEstencil> &consultas,
                                std::vector<Estencil> &saida) const
{
  for (size_t q = 0; q < consultas.size(); q++) {
    const ConsultaEstencil &c = consultas[q];
    saida[q].ids = c.ids;
    saida[q].pesos = mls_.pesos(c.alvo, c.suporte);
    saida[q].rhs = 0.0;   // sem fechamento de contorno ainda -- ver o cabecalho
  }
}

#ifdef T8FLOW_COM_OPENCL

void LoteEstenceis::calcula_gpu(const std::vector<ConsultaEstencil> &consultas,
                                std::vector<Estencil> &saida) const
{
  // O kernel exige suporte de tamanho UNIFORME no lote -- um item de trabalho por
  // consulta, com o mesmo `ksup`.  Lote heterogeneo cai na CPU inteiro, de
  // proposito: partir o lote por tamanho e' otimizacao, e otimizacao sem medida e'
  // onde o erro entra.
  const size_t ksup = consultas[0].suporte.size();
  for (const ConsultaEstencil &c : consultas)
    if (c.suporte.size() != ksup) { calcula_cpu(consultas, saida); return; }

  const size_t n = consultas.size();
  const bool d64 = gpu_->fp64;
  const size_t tam = d64 ? sizeof(double) : sizeof(float);

  std::vector<double> sup_d, alvo_d;
  std::vector<float> sup_f, alvo_f;
  if (d64) {
    sup_d.reserve(n * ksup * 2); alvo_d.reserve(n * 2);
    for (const ConsultaEstencil &c : consultas) {
      for (const Ponto2 &p : c.suporte) { sup_d.push_back(p[0]); sup_d.push_back(p[1]); }
      alvo_d.push_back(c.alvo[0]); alvo_d.push_back(c.alvo[1]);
    }
  } else {
    sup_f.reserve(n * ksup * 2); alvo_f.reserve(n * 2);
    for (const ConsultaEstencil &c : consultas) {
      for (const Ponto2 &p : c.suporte) {
        sup_f.push_back(static_cast<float>(p[0]));
        sup_f.push_back(static_cast<float>(p[1]));
      }
      alvo_f.push_back(static_cast<float>(c.alvo[0]));
      alvo_f.push_back(static_cast<float>(c.alvo[1]));
    }
  }

  cl_int e = CL_SUCCESS;
  void *psup = d64 ? static_cast<void *>(sup_d.data()) : static_cast<void *>(sup_f.data());
  void *palvo = d64 ? static_cast<void *>(alvo_d.data()) : static_cast<void *>(alvo_f.data());

  cl_mem bsup = clCreateBuffer(gpu_->ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                               n * ksup * 2 * tam, psup, &e);
  cl_mem balvo = clCreateBuffer(gpu_->ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                n * 2 * tam, palvo, &e);
  cl_mem bout = clCreateBuffer(gpu_->ctx, CL_MEM_WRITE_ONLY, n * ksup * tam, nullptr, &e);
  if (e != CL_SUCCESS) { calcula_cpu(consultas, saida); return; }

  const int ks = static_cast<int>(ksup);
  clSetKernelArg(gpu_->kernel, 0, sizeof(int), &ks);
  clSetKernelArg(gpu_->kernel, 1, sizeof(cl_mem), &bsup);
  clSetKernelArg(gpu_->kernel, 2, sizeof(cl_mem), &balvo);
  clSetKernelArg(gpu_->kernel, 3, sizeof(cl_mem), &bout);

  size_t global = n;
  e = clEnqueueNDRangeKernel(gpu_->fila, gpu_->kernel, 1, nullptr, &global, nullptr,
                             0, nullptr, nullptr);
  if (e == CL_SUCCESS) e = clFinish(gpu_->fila);

  if (e == CL_SUCCESS) {
    std::vector<double> out_d;
    std::vector<float> out_f;
    if (d64) { out_d.resize(n * ksup); e = clEnqueueReadBuffer(gpu_->fila, bout, CL_TRUE, 0, n*ksup*tam, out_d.data(), 0, nullptr, nullptr); }
    else     { out_f.resize(n * ksup); e = clEnqueueReadBuffer(gpu_->fila, bout, CL_TRUE, 0, n*ksup*tam, out_f.data(), 0, nullptr, nullptr); }

    if (e == CL_SUCCESS) {
      for (size_t q = 0; q < n; q++) {
        saida[q].ids = consultas[q].ids;
        saida[q].pesos.resize(ksup);
        for (size_t k = 0; k < ksup; k++)
          saida[q].pesos[k] = d64 ? out_d[q*ksup+k] : static_cast<double>(out_f[q*ksup+k]);
        saida[q].rhs = 0.0;
      }
    }
  }

  clReleaseMemObject(bsup); clReleaseMemObject(balvo); clReleaseMemObject(bout);

  // FALHA DA PLACA CAI NA CPU, mas nao em silencio quanto ao resultado: o que cai
  // e' o CALCULO, e o valor devolvido continua sendo o correto.  Esconder a falha
  // seria devolver lixo; esconder o fallback seria perder desempenho sem avisar --
  // e por isso o motivo vai para `porque_nao`, consultavel.
  if (e != CL_SUCCESS) {
    gpu_->porque_nao = "a execucao na placa falhou (codigo " + std::to_string(e) +
                       "); o lote foi calculado na CPU";
    calcula_cpu(consultas, saida);
  }
}

#else

void LoteEstenceis::calcula_gpu(const std::vector<ConsultaEstencil> &consultas,
                                std::vector<Estencil> &saida) const
{
  // Sem OpenCL em tempo de compilacao nao ha' como chegar aqui -- `usa_back_end`
  // lanca antes.  A definicao existe para o arquivo compilar igual nos dois casos.
  calcula_cpu(consultas, saida);
}

#endif

}  // namespace t8flow2d
