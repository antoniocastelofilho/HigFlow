// O LEITOR DE CONFIGURACAO, verificado contra CASOS REAIS DO HIGFLOW.
//
// O ORACULO CENTRAL DESTE ARQUIVO NAO E' UM ARQUIVO QUE EU ESCREVI: e' o caso
// `example2d_Newt` do HiGFlow, lido do lugar onde ele vive.  Um teste sobre
// YAML de brinquedo prova que o leitor le' o que eu imaginei do esquema; ler o
// caso de verdade prova que ele le' O ESQUEMA.
//
// E e' por isso que o caminho do HiGFlow entra aqui como opcional: se a arvore
// dele nao estiver ao lado, o caso e' PULADO COM AVISO -- nao silenciosamente
// aprovado, que seria o jeito de um teste verde nao provar nada.

#include "../testing.h"

#include <cstdio>
#include <stdexcept>
#include <string>

#include "configuracao.hxx"

using namespace t8flow;

namespace {

const char *HIGFLOW = "/home/castelo/HiGFlow-System/HigFlow/higflow";

bool existe(const std::string &c)
{
  std::FILE *f = std::fopen(c.c_str(), "r");
  if (f == nullptr) return false;
  std::fclose(f);
  return true;
}

std::string escreve(const std::string &nome, const std::string &texto)
{
  const std::string caminho = "/tmp/t8flow-cfg-" + nome;
  std::FILE *f = std::fopen(caminho.c_str(), "w");
  std::fputs(texto.c_str(), f);
  std::fclose(f);
  return caminho;
}

}  // namespace

int main(int argc, char **argv)
{
  (void) argc; (void) argv;

  t_case("le_o_caso_newtoniano_do_higflow_sem_traducao");
  {
    const std::string pre = std::string(HIGFLOW) +
                            "/example2d_Newt/input/example-2d.load";
    if (!existe(pre + ".par.contr.yaml")) {
      T_PULADO("a arvore do HiGFlow nao esta' em %s -- o caso REAL nao foi lido, "
             "e este arquivo so' verificou o esquema que eu mesmo escrevi", HIGFLOW);
    } else {
      Configuracao c;
      c.junta_caso(pre);

      // OS VALORES SAO OS DO ARQUIVO, e estao afirmados aqui para que uma
      // mudanca no caso do HiGFlow apareca como falha e nao como deriva.
      T_CHECK_MSG(c.inteiro("simulation_par.numsteps") == 101,
                  "numsteps devia ser 101, veio %ld",
                  c.inteiro("simulation_par.numsteps"));
      T_CHECK_MSG(std::abs(c.real("simulation_par.dt") - 0.001) < 1e-15,
                  "dt devia ser 0.001, veio %.17g", c.real("simulation_par.dt"));
      T_CHECK_MSG(std::abs(c.real("adimensional.Re") - 1.0) < 1e-15,
                  "Re devia ser 1.0");
      T_CHECK_MSG(c.texto("simulation_contr.projtype") == "incremental",
                  "projtype veio \"%s\"",
                  c.texto("simulation_contr.projtype").c_str());
      T_CHECK_MSG(c.texto("flowphase") == "singlephase",
                  "flowphase veio \"%s\"", c.texto("flowphase").c_str());
      T_CHECK_MSG(c.texto("singlephase.contr.flowtype") == "newtonian",
                  "flowtype veio \"%s\"",
                  c.texto("singlephase.contr.flowtype").c_str());

      // O ESQUEMA E' ANINHADO DE VERDADE, e o caminho pontilhado tem de chegar
      // la': tres niveis abaixo da raiz.
      T_CHECK_MSG(c.tem("singlephase.viscoelastic.model_giesekus.alpha"),
                  "o caminho de tres niveis nao chegou ao model_giesekus");

      // E O BOOLEANO DO YAML nao e' texto.
      T_CHECK_MSG(c.booleano("singlephase.contr.eoflow") == false,
                  "eoflow devia ser false");
    }
  }

  t_case("a_chave_que_falta_lanca_dizendo_o_bloco");
  {
    const std::string a = escreve("falta.par.contr.yaml",
                                  "simulation_par:\n  dt: 0.01\n  numsteps: 10\n");
    Configuracao c;
    c.junta(a);

    bool lancou = false;
    std::string msg;
    try { c.real("simulation_par.dts"); }
    catch (const std::runtime_error &e) { lancou = true; msg = e.what(); }
    T_CHECK_MSG(lancou, "chave obrigatoria ausente foi aceita");
    // A MENSAGEM TEM DE SERVIR: sem as vizinhas, quem le' nao sabe se errou o
    // nome ou o nivel.
    T_CHECK_MSG(msg.find("simulation_par.dt") != std::string::npos,
                "a mensagem nao lista as chaves do mesmo bloco: \"%s\"", msg.c_str());
  }

  t_case("o_nome_de_opcao_errado_nao_cai_no_padrao");
  {
    // O CASO QUE MOTIVA A FUNCAO.  Um `b` a mais em "cubista" nao casa com
    // nenhum esquema conhecido; um leitor permissivo roda o padrao e a corrida
    // sai com OUTRA discretizacao, plausivel e errada.
    const std::string a = escreve("typo.par.contr.yaml",
                                  "simulation_contr:\n"
                                  "  secondconvecdiscrtype: cubbista\n");
    Configuracao c;
    c.junta(a);
    const std::vector<std::string> validos = {"modified_coefficient_upwind",
                                              "cubista", "quick"};
    bool lancou = false;
    std::string msg;
    try { c.opcao("simulation_contr.secondconvecdiscrtype", validos, 1); }
    catch (const std::runtime_error &e) { lancou = true; msg = e.what(); }
    T_CHECK_MSG(lancou, "\"cubbista\" foi aceito e caiu no padrao");
    T_CHECK_MSG(msg.find("cubista") != std::string::npos,
                "a mensagem nao traz a lista de validos: \"%s\"", msg.c_str());

    // CONTRASTE: o nome certo passa, e devolve o indice certo.  Sem este par,
    // "lanca" nao distingue um leitor correto de um que recusa tudo.
    const std::string b = escreve("certo.par.contr.yaml",
                                  "simulation_contr:\n"
                                  "  secondconvecdiscrtype: quick\n");
    Configuracao d;
    d.junta(b);
    T_CHECK_MSG(d.opcao("simulation_contr.secondconvecdiscrtype", validos) == 2,
                "\"quick\" devia dar indice 2");

    // E A CHAVE AUSENTE usa o padrao de quem chama, que esta' na chamada.
    Configuracao vazia;
    T_CHECK_MSG(vazia.opcao("simulation_contr.secondconvecdiscrtype", validos, 1) == 1,
                "o padrao explicito nao foi usado");
  }

  t_case("a_chave_nao_lida_e_denunciada");
  {
    // O MODO DE FALHA QUE ISTO FECHA: `numstpes` com as letras trocadas nao e'
    // lido por ninguem, a corrida usa o padrao do codigo, e nada acusa.
    const std::string a = escreve("sobra.par.contr.yaml",
                                  "simulation_par:\n"
                                  "  dt: 0.01\n"
                                  "  numstpes: 10\n");
    Configuracao c;
    c.junta(a);
    c.real("simulation_par.dt");

    const std::vector<std::string> fora = c.nao_lidas();
    T_CHECK_MSG(fora.size() == 1, "deviam sobrar 1 chave, sobraram %zu", fora.size());
    T_CHECK_MSG(fora.size() == 1 && fora[0] == "simulation_par.numstpes",
                "a chave denunciada devia ser numstpes");

    // CONTRASTE: ler tudo nao deixa sobra.
    Configuracao d;
    d.junta(a);
    d.real("simulation_par.dt");
    d.inteiro("simulation_par.numstpes");
    T_CHECK_MSG(d.nao_lidas().empty(), "nada devia sobrar quando tudo foi lido");
  }

  t_case("o_tipo_e_cobrado");
  {
    const std::string a = escreve("tipo.par.contr.yaml",
                                  "simulation_par:\n"
                                  "  dt: muito pouco\n"
                                  "  numsteps: 10.5\n"
                                  "  gravar: talvez\n");
    Configuracao c;
    c.junta(a);

    // `strtod` sozinho devolveria 0 para "muito pouco", e zero e' um dt
    // plausivel de digitar -- a corrida rodaria parada em vez de recusar.
    bool r = false, i = false, b = false;
    try { c.real("simulation_par.dt"); } catch (const std::runtime_error &) { r = true; }
    try { c.inteiro("simulation_par.numsteps"); } catch (const std::runtime_error &) { i = true; }
    try { c.booleano("simulation_par.gravar"); } catch (const std::runtime_error &) { b = true; }
    T_CHECK_MSG(r, "\"muito pouco\" foi aceito como real");
    T_CHECK_MSG(i, "\"10.5\" foi aceito como inteiro");
    T_CHECK_MSG(b, "\"talvez\" foi aceito como booleano");
  }

  t_case("a_colisao_entre_arquivos_do_caso_e_erro");
  {
    // DOIS ARQUIVOS DO MESMO CASO declarando a mesma chave com valores
    // diferentes e' defeito do caso.  O ultimo a ser lido vencer faz o
    // resultado depender da ordem de leitura, que ninguem escreveu.
    const std::string a = escreve("col1.par.contr.yaml", "simulation_par:\n  dt: 0.01\n");
    Configuracao c;
    c.junta(a);
    const std::string b = escreve("col2.par.contr.yaml", "simulation_par:\n  dt: 0.02\n");
    bool lancou = false;
    try { c.junta(b); } catch (const std::runtime_error &) { lancou = true; }
    T_CHECK_MSG(lancou, "a colisao com valor diferente foi aceita em silencio");

    // CONTRASTE: o MESMO valor nos dois nao e' colisao -- e' redundancia, que
    // acontece em caso remendado e nao muda resultado nenhum.
    Configuracao d;
    d.junta(a);
    bool ok = true;
    try { d.junta(escreve("col3.par.contr.yaml", "simulation_par:\n  dt: 0.01\n")); }
    catch (const std::runtime_error &) { ok = false; }
    T_CHECK_MSG(ok, "o mesmo valor repetido foi tratado como colisao");
  }

  return t_end();
}
