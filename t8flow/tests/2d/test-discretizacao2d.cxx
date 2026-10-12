// O modulo de escolha entre DIFERENCAS FINITAS e VOLUMES FINITOS.
//
// O QUE ESTE ARREIO PROTEGE: a coerencia entre a familia escolhida e as quatro
// consequencias dela.  Separa-las e' um erro que NAO FALHA -- um operador de
// volumes finitos com fonte pontual continua convergindo, em ordem reduzida e
// com a constante errada --, e por isso nenhum teste de convergencia o pega.
//
// Cobre as clausulas W1, W2 e W6 do contrato.

#include "../testing.h"

#include <stdexcept>
#include <string>

#include "arranque.hxx"
#include "discretizacao2d.hxx"

using namespace t8flow2d;

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  t_case("a_familia_sai_da_variante_e_o_nome_traz_o_prefixo");
  {
    // A TABELA INTEIRA, de proposito: uma variante nova que entre no `enum` sem
    // entrar aqui passa despercebida, e o comentario em `familia_de` confia no
    // aviso do compilador -- que so' cobre o `switch`, nao a semantica.
    const struct { Variante v; Familia f; const char *nome; } tab[] = {
      {Variante::FluxoBalanco,        Familia::VolumesFinitos,    "vf-fluxo2"},
      {Variante::MlsGauss,            Familia::VolumesFinitos,    "vf-mlsgauss3"},
      {Variante::CompactoMehrstellen, Familia::DiferencasFinitas, "df-mehrstellen4"},
      {Variante::CompactoLele,        Familia::DiferencasFinitas, "df-lele4"},
      {Variante::MlsPontual,          Familia::DiferencasFinitas, "df-mlspontual3"},
    };
    for (const auto &t : tab) {
      const Discretizacao d(t.v, 3);
      T_CHECK_MSG(d.familia() == t.f, "familia errada para \"%s\"", t.nome);
      T_CHECK_MSG(d.nome() == t.nome, "nome \"%s\", esperado \"%s\"",
                  d.nome().c_str(), t.nome);

      // O PREFIXO DO NOME TEM DE CASAR COM A FAMILIA.  Quem escolhe na linha de
      // comando ve' o nome, nao o `enum` -- se o prefixo mentir, a escolha mente.
      const std::string pre = d.nome().substr(0, 3);
      T_CHECK_MSG(pre == (d.familia() == Familia::VolumesFinitos ? "vf-" : "df-"),
                  "o prefixo de \"%s\" nao corresponde a' familia", t.nome);

      // E a ida e volta pelo nome preserva a escolha.
      const Discretizacao v = Discretizacao::do_nome(t.nome, 3);
      T_CHECK_MSG(v.variante() == t.v && v.familia() == t.f,
                  "ida e volta por \"%s\" mudou a escolha", t.nome);
    }
    // A ORDEM E' PARTE DA ESCOLHA, e o nome tem de carrega-la: "volumes finitos
    // de quarta ordem" precisa ser nomeavel, nao um parametro separado que
    // alguem tenha de descobrir que existe.
    for (int o : {2, 3, 4}) {
      const Discretizacao d(Variante::MlsGauss, o);
      const std::string esperado = "vf-mlsgauss" + std::to_string(o);
      T_CHECK_MSG(d.nome() == esperado, "nome \"%s\", esperado \"%s\"",
                  d.nome().c_str(), esperado.c_str());
      T_CHECK_MSG(d.ordem() == o, "ordem nominal %d, esperada %d", d.ordem(), o);
      const Discretizacao v = Discretizacao::do_nome(esperado);
      T_CHECK_MSG(v.ordem_mls() == o,
                  "\"%s\" voltou com ordem %d", esperado.c_str(), v.ordem_mls());
    }

    // E A FORMA CURTA CONTINUA VALENDO: `vf-mlsgauss` sem sufixo toma a ordem de
    // quem chama.  Tira-la quebraria quem ja' a usa sem ganhar nada.
    T_CHECK_MSG(Discretizacao::do_nome("vf-mlsgauss", 4).ordem_mls() == 4,
                "a forma curta nao respeitou a ordem informada");

    // O `2` DE `vf-fluxo2` NAO E' SUFIXO DE ORDEM -- e' parte do nome.  Sem esta
    // distincao, "vf-fluxo2" viraria "vf-fluxo" de ordem 2 e nao resolveria.
    T_CHECK_MSG(Discretizacao::do_nome("vf-fluxo2").variante() == Variante::FluxoBalanco,
                "o digito de \"vf-fluxo2\" foi lido como sufixo de ordem");

    printf("  (W1: 5 variantes mais as ordens de FV, familia e nome conferidos)\n");
  }

  t_case("as_consequencias_da_familia_andam_juntas");
  {
    // VOLUMES FINITOS: equacao integrada, fonte integrada, residuo pela area.
    //
    // NAO "incognita media": os dois eixos sao independentes, e nenhuma variante
    // deste sistema usa incognita-media.  Medido: `vf-mlsgauss` de ordem 3 e'
    // exato para grau 3 comparado com o VALOR NO CENTRO, o que nao aconteceria
    // se a incognita fosse a media -- a diferenca entre a media e o valor
    // central de um cubico nao e' zero.
    for (Variante v : {Variante::FluxoBalanco, Variante::MlsGauss}) {
      const Discretizacao d(v, 3);
      T_CHECK_MSG(d.equacao_integrada(), "%s devia integrar a equacao",
                  d.nome().c_str());
      T_CHECK_MSG(!d.incognita_e_media(),
                  "%s foi declarado com incognita MEDIA; nenhuma variante deste "
                  "sistema usa, e o oraculo compara com o valor no centro",
                  d.nome().c_str());
      T_CHECK_MSG(d.fonte_integrada(), "%s devia integrar a fonte", d.nome().c_str());
      T_NEAR(d.normalizacao_do_residuo(0.25), 0.25, 1e-15,
             "normalizacao do residuo em volumes finitos");
    }

    // DIFERENCAS FINITAS: incognita e' valor, fonte pontual, residuo sem dividir.
    for (Variante v : {Variante::CompactoMehrstellen, Variante::CompactoLele,
                       Variante::MlsPontual}) {
      const Discretizacao d(v, 3);
      T_CHECK_MSG(!d.equacao_integrada(), "%s nao devia integrar a equacao",
                  d.nome().c_str());
      T_CHECK_MSG(!d.incognita_e_media(), "%s devia ter incognita pontual",
                  d.nome().c_str());
      T_CHECK_MSG(!d.fonte_integrada(), "%s nao devia integrar a fonte",
                  d.nome().c_str());
      // DIVIDIR PELA AREA EM DF daria um numero que CRESCE com o refino sem
      // significado, e esse foi um engano real deste projeto.
      T_NEAR(d.normalizacao_do_residuo(0.25), 1.0, 1e-15,
             "normalizacao do residuo em diferencas finitas");
    }
    printf("  (W2: as quatro consequencias conferidas nas 5 variantes)\n");
  }

  t_case("a_instabilidade_do_pontual_em_malha_adaptada_e_declarada");
  {
    // O QUE SE PROIBE AQUI E' A RESPOSTA OTIMISTA.  Medido: o laplaciano pontual
    // com ordem 3 em malha adaptada nao converge em 20 mil iteracoes, e o erro
    // vai a 1,3e+06.  Declarar que trata nao basta -- tem de declarar que nao e'
    // estavel, que e' outra pergunta.
    const Discretizacao pont(Variante::MlsPontual, 3);
    T_CHECK_MSG(pont.trata_malha_adaptada(),
                "o laplaciano pontual TRATA malha adaptada -- o operador esta' "
                "definido nela, ao contrario dos compactos");
    T_CHECK_MSG(!pont.estavel_em_malha_adaptada(),
                "o laplaciano pontual foi declarado ESTAVEL em malha adaptada; "
                "medido, ele nao converge com ordem 3");

    // O CONTRASTE: o balanco de fluxos garante a estrutura por construcao, e as
    // duas variantes de volumes finitos sao estaveis.  Sem este par, "nao e'
    // estavel" nao distingue uma propriedade do metodo de um defeito geral.
    for (Variante v : {Variante::FluxoBalanco, Variante::MlsGauss}) {
      const Discretizacao d(v, 3);
      T_CHECK_MSG(d.trata_malha_adaptada() && d.estavel_em_malha_adaptada(),
                  "%s devia tratar E ser estavel em malha adaptada", d.nome().c_str());
    }

    // E as compactas recusam antes: nao tratam, logo a pergunta da estabilidade
    // nem se poe.
    for (Variante v : {Variante::CompactoMehrstellen, Variante::CompactoLele}) {
      const Discretizacao d(v);
      T_CHECK_MSG(!d.trata_malha_adaptada(),
                  "%s devia RECUSAR malha adaptada", d.nome().c_str());
    }
    printf("  (W6: pontual trata mas nao e' estavel; fluxo trata e e')\n");
  }

  t_case("volumes_finitos_acima_de_segunda_ordem_existe_e_e_nomeavel");
  {
    // A PERGUNTA QUE ESTE CASO RESPONDE: "tem como discretizar por FV com ordem
    // maior que 2?".  Tem -- `vf-mlsgauss3` e `vf-mlsgauss4` --, e o caso afirma
    // as tres coisas que tornam isso uma resposta e nao uma promessa: a opcao e'
    // nomeavel, e' de volumes finitos, e trata malha adaptada.
    for (int o : {3, 4}) {
      const Discretizacao d = Discretizacao::do_nome("vf-mlsgauss" + std::to_string(o));
      T_CHECK_MSG(d.familia() == Familia::VolumesFinitos,
                  "vf-mlsgauss%d nao e' de volumes finitos", o);
      T_CHECK_MSG(d.equacao_integrada(),
                  "vf-mlsgauss%d nao integra a equacao -- sem isso nao e' "
                  "conservativo, que e' a razao de escolher volumes finitos", o);
      T_CHECK_MSG(d.ordem() > 2, "vf-mlsgauss%d tem ordem nominal %d", o, d.ordem());
      T_CHECK_MSG(d.trata_malha_adaptada() && d.estavel_em_malha_adaptada(),
                  "vf-mlsgauss%d nao serve em malha adaptada", o);
      T_CHECK_MSG(d.grau_exato(true) == o && d.grau_exato(false) == o,
                  "vf-mlsgauss%d: grau %d/%d, esperado %d nos dois tipos de malha",
                  o, d.grau_exato(false), d.grau_exato(true), o);
    }

    // O CONTRASTE, que e' o que da' sentido ao caso: a outra variante de volumes
    // finitos NAO passa de segunda ordem, e perde exatidao em malha adaptada.
    const Discretizacao f(Variante::FluxoBalanco);
    T_CHECK_MSG(f.ordem() == 2 && f.grau_exato(true) == 1,
                "o fluxo de segunda ordem deveria ter ordem 2 e grau 1 em malha "
                "adaptada; deu %d e %d", f.ordem(), f.grau_exato(true));
    printf("  (FV acima de 2a ordem: vf-mlsgauss3 e vf-mlsgauss4, exatos ate' o "
           "grau da ordem nas duas malhas)\n");
  }

  t_case("o_grau_de_exatidao_e_por_malha_e_vem_de_medida");
  {
    // O GRAU MUDA COM O TIPO DE MALHA no esquema de fluxo, e nao muda nas
    // variantes de MLS -- e e' justamente essa diferenca que faz as de MLS
    // existirem.  Os numeros vem do relatorio, nao de deducao.
    const Discretizacao f(Variante::FluxoBalanco);
    T_CHECK_MSG(f.grau_exato(false) == 2 && f.grau_exato(true) == 1,
                "fluxo: grau %d na uniforme e %d na adaptada (esperado 2 e 1)",
                f.grau_exato(false), f.grau_exato(true));

    for (int o : {2, 3, 4}) {
      const Discretizacao g(Variante::MlsGauss, o);
      T_CHECK_MSG(g.grau_exato(false) == o && g.grau_exato(true) == o,
                  "mlsgauss ordem %d: grau %d/%d, esperado %d nos dois",
                  o, g.grau_exato(false), g.grau_exato(true), o);
    }

    const Discretizacao c(Variante::CompactoLele);
    T_CHECK_MSG(c.grau_exato(false) == 3,
                "compacto: grau %d, esperado 3 -- o limite e' o fechamento de "
                "parede, nao o estencil interior", c.grau_exato(false));
    printf("  (grau: fluxo 2/1, mls = ordem, compacto 3)\n");
  }

  t_case("escolha_malformada_e_recusada_com_a_lista");
  {
    bool lancou = false;
    std::string msg;
    try {
      Discretizacao::do_nome("volumes-finitos");
    } catch (const std::invalid_argument &e) {
      lancou = true;
      msg = e.what();
    }
    T_CHECK_MSG(lancou, "nome invalido passou sem reclamar");
    // A LISTA VAI NA MENSAGEM: um erro que diz apenas "invalido" obriga quem o
    // recebe a ler o codigo.
    T_CHECK_MSG(msg.find("vf-fluxo2") != std::string::npos &&
                    msg.find("df-mlspontual") != std::string::npos,
                "a mensagem nao traz a lista de validos: \"%s\"", msg.c_str());

    // ORDEM DE MLS FORA DA FAIXA tambem e' recusada -- 4 e' o teto da base.
    bool lancou_ordem = false;
    try {
      Discretizacao(Variante::MlsGauss, 5);
    } catch (const std::invalid_argument &) { lancou_ordem = true; }
    T_CHECK_MSG(lancou_ordem, "ordem 5 do MLS passou; o teto da base propria e' 4");

    // E a ordem NAO e' lida nas variantes que nao usam MLS: vale zero, e isso e'
    // informacao -- um valor diferente de zero ali significaria que alguem
    // acredita estar configurando algo que o esquema nao le'.
    T_CHECK_MSG(Discretizacao(Variante::CompactoLele, 3).ordem_mls() == 0,
                "a variante compacta reportou ordem de MLS diferente de zero");
    printf("  (recusas: nome, ordem fora da faixa, ordem em variante que nao usa)\n");
  }

  t_case("o_arranjo_e_deslocado_e_cada_operador_diz_de_onde_le_e_para_onde_escreve");
  {
    // O ARRANJO DESLOCADO (MAC) NAS DUAS FAMILIAS -- decisao do usuario.  Ela so'
    // significa algo se cada operador disser as suas duas localizacoes, e e'
    // isso que este caso afirma.
    for (Variante v : {Variante::FluxoBalanco, Variante::MlsGauss,
                       Variante::CompactoLele, Variante::MlsPontual}) {
      T_CHECK_MSG(Discretizacao(v, 3).arranjo_deslocado(),
                  "%s nao declarou arranjo deslocado", Discretizacao(v, 3).nome().c_str());
    }

    for (int d = 0; d < 2; d++) {
      const Local fac = (d == 0) ? Local::FacetaX : Local::FacetaY;

      // O GRADIENTE DE PRESSAO E O DIVERGENTE SAO DUAIS: um vai do centro para a
      // faceta, o outro da faceta para o centro.  E' dessa dualidade que vem a
      // conservacao de massa na projecao, e troca-la seria trocar o esquema.
      T_CHECK_MSG(origem_de(Operador::GradientePressao, d) == Local::Centro &&
                      destino_de(Operador::GradientePressao, d) == fac,
                  "grad p na direcao %d nao vai do centro para a faceta", d);
      T_CHECK_MSG(origem_de(Operador::DivergenteVelocidade, d) == fac &&
                      destino_de(Operador::DivergenteVelocidade, d) == Local::Centro,
                  "div u na direcao %d nao vai da faceta para o centro", d);

      // O VISCOSO E A ADVECCAO FICAM NA FACETA, entrada e saida.
      for (Operador op : {Operador::LaplacianoVelocidade, Operador::Adveccao})
        T_CHECK_MSG(origem_de(op, d) == fac && destino_de(op, d) == fac,
                    "operador de faceta na direcao %d saiu do lugar", d);
    }

    // SO' A ADVECCAO INTERPOLA, e e' o unico preco do arranjo deslocado: para
    // avancar a componente `x` e' preciso da componente `y` sobre a faceta de
    // `x`, e elas nao moram juntas.  Afirmar que os outros tres NAO interpolam e'
    // a outra metade -- e' isso que o MAC compra.
    T_CHECK_MSG(exige_interpolacao(Operador::Adveccao),
                "a adveccao deveria exigir interpolacao");
    for (Operador op : {Operador::GradientePressao, Operador::DivergenteVelocidade,
                        Operador::LaplacianoVelocidade})
      T_CHECK_MSG(!exige_interpolacao(op),
                  "um operador que nao deveria interpolar declarou que interpola");

    // DIRECAO INVALIDA E' RECUSADA: em 2D so' ha' 0 e 1, e aceitar 2 em silencio
    // devolveria a faceta de `y` para uma direcao que nao existe.
    bool lancou = false;
    try { origem_de(Operador::GradientePressao, 2); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "direcao 2 passou num codigo de duas dimensoes");

    printf("  (arranjo MAC: grad centro->faceta, div faceta->centro, "
           "viscoso e adveccao na faceta; so' a adveccao interpola)\n");
  }

  return t_end();
}
