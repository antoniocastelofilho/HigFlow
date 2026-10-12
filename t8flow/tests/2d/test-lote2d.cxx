// O LOTE DE ESTENCEIS, com os dois back-ends.
//
// A CPU E' AFIRMADA SEMPRE; a GPU so' quando ha' uma.  E isso nao e' concessao:
// a versao de CPU e' a implementacao de primeira classe, e quem roda sem placa
// tem de ter a suite INTEIRA verde.  Por isso nenhum caso aqui depende de GPU
// para passar -- os casos de placa se declaram SALTADOS quando nao ha' placa, e
// um salto e' visivel na saida, nao silencioso.
//
// ORACULOS:
//
//   reproducao com contraste   o mesmo da suite: um lote de ordem 1 reproduz
//                              campo linear exatamente e erra no quadratico.  Vale
//                              para QUALQUER back-end, e e' o que define correto;
//   soma dos pesos = 1         reproducao da constante, independente do campo --
//                              pega erro de escala que a reproducao de um
//                              polinomio especifico pode absorver;
//   GPU contra CPU             ponto a ponto, nos PESOS e nao no valor
//                              interpolado: o valor pode bater por cancelamento
//                              com pesos errados, os pesos nao;
//   recusa igual nos dois      suporte invalido e' recusado antes de escolher o
//                              back-end, senao "funciona na CPU e quebra na placa"
//                              passa a ser possivel.
#include "../testing.h"

#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "lote2d.hxx"

using t8flow2d::BackEnd;
using t8flow2d::ConsultaEstencil;
using t8flow2d::Estencil;
using t8flow2d::LoteEstenceis;
using t8flow2d::Mls;
using t8flow2d::Ponto2;

static double grau1(const Ponto2 &p) { return 1.0 + 2.0 * p[0] - 3.0 * p[1]; }
static double grau2(const Ponto2 &p) { return grau1(p) + 4.0 * p[0] * p[0]; }

//! Lote deterministico: alvos assimetricos de proposito -- no centro de simetria
//! do suporte os termos impares se cancelam e um ajuste de ordem 1 reproduz
//! quadratica por acidente.  Medido na suite, vale aqui igual.
static std::vector<ConsultaEstencil> monta_lote(int ncons, int ksup, double h)
{
  std::vector<ConsultaEstencil> v;
  v.reserve(static_cast<size_t>(ncons));
  for (int q = 0; q < ncons; q++) {
    ConsultaEstencil c;
    c.alvo = {0.37 + 0.001 * (q % 211), 0.44 + 0.001 * (q % 307)};
    for (int k = 0; k < ksup; k++) {
      c.suporte.push_back({c.alvo[0] + ((k % 4) - 1.3) * h,
                           c.alvo[1] + ((k / 4) - 1.7) * h});
      c.ids.push_back(1000 + q * ksup + k);
    }
    v.push_back(c);
  }
  return v;
}

static double interpola(const Estencil &e, const ConsultaEstencil &c,
                        double (*f)(const Ponto2 &))
{
  double v = 0.0;
  for (size_t k = 0; k < e.pesos.size(); k++) v += e.pesos[k] * f(c.suporte[k]);
  return v;
}

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  const int NCONS = 2000, KSUP = 16;
  const std::vector<ConsultaEstencil> lote = monta_lote(NCONS, KSUP, 0.1);

  LoteEstenceis L(1, KSUP);
  std::vector<Estencil> cpu;

  // O PADRAO SEGUE QUEM COMPILOU.  Com `make GPU=sim` a GPU e' o padrao e o lote
  // ja' comeca nela; nos outros modos comeca na CPU.  O caso afirma a COERENCIA
  // entre a flag e o comportamento, que e' o que se pode afirmar de dentro de um
  // unico build -- e imprime o estado, para a saida dizer qual build foi testado.
  t_case("o_back_end_padrao_segue_a_flag_de_compilacao");
  {
    const bool pedida = LoteEstenceis::gpu_pedida_na_compilacao();
    printf("  (compilado com GPU padrao: %s; GPU disponivel: %s)\n",
           pedida ? "sim" : "nao", L.gpu_disponivel() ? "sim" : "nao");
    if (pedida && L.gpu_disponivel()) {
      T_CHECK_MSG(L.back_end() == BackEnd::Gpu,
                  "compilado com GPU=sim e com placa, mas o padrao nao e' GPU");
    } else {
      T_CHECK_MSG(L.back_end() == BackEnd::Cpu,
                  "sem GPU=sim (ou sem placa), o padrao tem de ser CPU");
    }
  }

  t_case("o_lote_de_cpu_devolve_um_estencil_por_consulta_com_ids_preservados");
  {
    // FORCA A CPU, e isto nao e' redundante: com `make GPU=sim` o padrao e' a
    // placa, e sem esta linha o vetor chamado `cpu` seria calculado NA GPU -- a
    // comparacao adiante compararia GPU contra GPU e daria diferenca exatamente
    // zero, que e' como o oraculo deixa de existir sem nada reclamar.  Foi o que
    // aconteceu: 0,000e+00 de diferenca denunciou.
    L.usa_back_end(BackEnd::Cpu);
    L.calcula(lote, cpu);
    T_CHECK_MSG(cpu.size() == lote.size(), "%zu estenceis para %zu consultas",
                cpu.size(), lote.size());
    bool ok = true;
    for (size_t q = 0; q < cpu.size() && ok; q++) {
      if (cpu[q].ids != lote[q].ids || cpu[q].pesos.size() != lote[q].suporte.size()) {
        T_CHECK_MSG(false, "consulta %zu: %zu ids e %zu pesos para %zu pontos",
                    q, cpu[q].ids.size(), cpu[q].pesos.size(), lote[q].suporte.size());
        ok = false;
      }
    }
  }

  t_case("cpu_reproduz_o_grau_da_ordem");
  double pior1 = 0.0, pior2 = 0.0;
  {
    for (size_t q = 0; q < cpu.size(); q++) {
      const double e1 = std::fabs(interpola(cpu[q], lote[q], grau1) - grau1(lote[q].alvo));
      if (e1 > pior1) pior1 = e1;
      const double e2 = std::fabs(interpola(cpu[q], lote[q], grau2) - grau2(lote[q].alvo));
      if (e2 > pior2) pior2 = e2;
    }
    T_NEAR(pior1, 0.0, 1e-12, "pior erro no campo de grau 1 (ajuste de ordem 1)");
  }

  t_case("cpu_erra_mais_no_grau_acima_da_ordem");
  T_CONTRASTE(pior1 > 1e-16 ? pior1 : 1e-16, pior2, 1e3,
              "grau 1 contra grau 2 no lote de ordem 1");

  t_case("os_pesos_somam_um");
  {
    double pior = 0.0;
    for (const Estencil &e : cpu) {
      double s = 0.0;
      for (double w : e.pesos) s += w;
      const double d = std::fabs(s - 1.0);
      if (d > pior) pior = d;
    }
    T_NEAR(pior, 0.0, 1e-12, "pior desvio da soma dos pesos");
  }

  // AS RECUSAS, antes de escolher o back-end.
  t_case("suporte_insuficiente_e_recusado_identificando_a_consulta");
  {
    std::vector<ConsultaEstencil> ruim = monta_lote(3, KSUP, 0.1);
    ruim[1].suporte.resize(2);
    ruim[1].ids.resize(2);
    bool lancou = false;
    try { std::vector<Estencil> s; L.calcula(ruim, s); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "consulta com 2 pontos passou sem reclamar");
  }

  t_case("ids_em_numero_diferente_do_suporte_e_recusado");
  {
    std::vector<ConsultaEstencil> ruim = monta_lote(2, KSUP, 0.1);
    ruim[0].ids.pop_back();
    bool lancou = false;
    try { std::vector<Estencil> s; L.calcula(ruim, s); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "ids em numero diferente do suporte passou sem reclamar");
  }

  // --- a GPU, se houver -----------------------------------------------------
  //
  // SEM PLACA NAO HA' FALHA, ha' SALTO -- e o salto aparece na saida.  Quem roda
  // sem OpenCL tem a suite inteira verde, que e' o ponto de ter versao de CPU.
  const bool tem_gpu = L.gpu_disponivel();

  t_case("sem_gpu_o_lote_continua_funcionando_e_diz_por_que");
  if (!tem_gpu) {
    const std::string d = L.diagnostico_gpu();
    T_CHECK_MSG(!d.empty(), "sem GPU e sem diagnostico -- o motivo tem de ser dizivel");
    printf("  (sem GPU: %s)\n", d.c_str());
    // E pedir GPU explicitamente tem de LANCAR, nao cair na CPU calado.
    bool lancou = false;
    try { L.usa_back_end(BackEnd::Gpu); }
    catch (const std::runtime_error &) { lancou = true; }
    T_CHECK_MSG(lancou, "pedir GPU sem GPU nao reclamou");
  } else {
    printf("  (GPU presente: %s)\n", L.diagnostico_gpu().c_str());
  }

  t_case("gpu_bate_com_cpu_nos_pesos");
  if (!tem_gpu) {
    printf("  SALTADO: nenhuma GPU utilizavel -- %s\n", L.diagnostico_gpu().c_str());
  } else {
    L.usa_back_end(BackEnd::Gpu);
    T_CHECK_MSG(L.back_end() == BackEnd::Gpu, "o back-end nao trocou para GPU");
    std::vector<Estencil> gpu;
    L.calcula(lote, gpu);
    T_CHECK_MSG(gpu.size() == cpu.size(), "%zu estenceis na GPU para %zu na CPU",
                gpu.size(), cpu.size());

    // A COMPARACAO E' NOS PESOS, nao no valor interpolado: valor pode bater por
    // cancelamento com pesos errados.
    //
    // E DIFERENCA EXATAMENTE ZERO E' SUSPEITA, nao boa noticia: significa que os
    // dois lados rodaram o MESMO codigo.  O caso recusa esse resultado, porque ele
    // indica referencia calculada na placa -- o modo como este teste ja' perdeu o
    // oraculo uma vez.
    double pior = 0.0;
    size_t pior_q = 0;
    for (size_t q = 0; q < gpu.size() && q < cpu.size(); q++)
      for (size_t k = 0; k < gpu[q].pesos.size() && k < cpu[q].pesos.size(); k++) {
        const double d = std::fabs(gpu[q].pesos[k] - cpu[q].pesos[k]);
        if (d > pior) { pior = d; pior_q = q; }
      }
    // A tolerancia depende da PRECISAO DO DISPOSITIVO, e e' por isso que o
    // diagnostico diz se ha' fp64: sem ela o kernel roda em float, e exigir 1e-12
    // de um dispositivo float seria reprovar a placa por uma escolha do driver.
    const bool float_so = L.diagnostico_gpu().find("SEM fp64") != std::string::npos;
    const double tol = float_so ? 1e-4 : 1e-10;
    T_CHECK_MSG(pior <= tol,
                "pior diferenca de peso %.3g na consulta %zu (tolerancia %.0e, "
                "dispositivo %s)", pior, pior_q, tol,
                float_so ? "em FLOAT" : "em double");
    T_CHECK_MSG(pior > 0.0,
                "diferenca EXATAMENTE zero entre GPU e CPU -- os dois lados rodaram "
                "o mesmo codigo, logo nao houve comparacao: a referencia nao foi "
                "calculada na CPU (o padrao e' GPU quando compilado com GPU=sim)");
    printf("  (GPU contra CPU: pior diferenca de peso %.3e, tolerancia %.0e)\n", pior, tol);

    // E a GPU tem de passar no MESMO oraculo da CPU -- bater com a CPU nao basta
    // se as duas estiverem erradas juntas.
    double g1 = 0.0, g2 = 0.0;
    for (size_t q = 0; q < gpu.size(); q++) {
      const double e1 = std::fabs(interpola(gpu[q], lote[q], grau1) - grau1(lote[q].alvo));
      if (e1 > g1) g1 = e1;
      const double e2 = std::fabs(interpola(gpu[q], lote[q], grau2) - grau2(lote[q].alvo));
      if (e2 > g2) g2 = e2;
    }
    T_CHECK_MSG(g1 <= tol, "a GPU nao reproduziu o grau 1: pior erro %.3g", g1);
    T_CONTRASTE(g1 > 1e-16 ? g1 : 1e-16, g2, 1e2,
                "grau 1 contra grau 2 na GPU");
    L.usa_back_end(BackEnd::Cpu);
  }

  t_case("lote_heterogeneo_cai_na_cpu_com_o_resultado_certo");
  if (!tem_gpu) {
    printf("  SALTADO: nenhuma GPU utilizavel\n");
  } else {
    // O kernel exige suporte uniforme; lote misturado cai na CPU INTEIRO, e o
    // resultado tem de continuar correto -- o que cai e' o calculo, nao o valor.
    std::vector<ConsultaEstencil> misto = monta_lote(10, KSUP, 0.1);
    misto[4].suporte.resize(KSUP - 4);
    misto[4].ids.resize(KSUP - 4);
    L.usa_back_end(BackEnd::Gpu);
    std::vector<Estencil> r;
    L.calcula(misto, r);
    double pior = 0.0;
    for (size_t q = 0; q < r.size(); q++) {
      const double e1 = std::fabs(interpola(r[q], misto[q], grau1) - grau1(misto[q].alvo));
      if (e1 > pior) pior = e1;
    }
    T_NEAR(pior, 0.0, 1e-12, "pior erro no lote heterogeneo (deve ter caido na CPU)");
    L.usa_back_end(BackEnd::Cpu);
  }

  t_case("ordem_que_o_kernel_nao_cobre_nao_vai_para_a_gpu_em_silencio");
  {
    // Ordem 3 em 2D sao 10 monomios; o kernel cobre 3 e 6.  A GPU tem de se
    // declarar INDISPONIVEL, nao calcular errado.
    LoteEstenceis L3(3, 40);
    T_CHECK_MSG(!L3.gpu_disponivel(),
                "ordem 3 apareceu como suportada pela GPU; diagnostico: %s",
                L3.diagnostico_gpu().c_str());
    // e a CPU atende normalmente
    std::vector<ConsultaEstencil> l3 = monta_lote(50, 36, 0.1);
    std::vector<Estencil> r;
    L3.calcula(l3, r);
    T_CHECK_MSG(r.size() == l3.size(), "a CPU nao atendeu a ordem 3");
  }

  // A SAIDA DIZ SE A GPU FOI EXERCITADA, e isso e' necessario: os casos de placa
  // PASSAM tanto comparando contra a CPU quanto saltando por nao haver placa, e o
  // verde sozinho nao distingue as duas coisas.  O driver colhe esta linha e a
  // reporta, para "contrato 100%" nunca significar "a GPU nao foi testada".
  printf("gpu-exercitada %s\n", tem_gpu ? "sim" : "nao");
  printf("gpu-padrao-de-compilacao %s\n",
         LoteEstenceis::gpu_pedida_na_compilacao() ? "sim" : "nao");
  return t_end();
}
