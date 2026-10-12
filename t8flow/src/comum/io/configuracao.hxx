// A CONFIGURACAO LIDA DE YAML, no mesmo esquema de arquivos do HiGFlow.
//
// POR QUE O MESMO ESQUEMA e nao um proprio, mais limpo: o objetivo declarado e'
// comparar os dois sistemas no mesmo caso, e um esquema proprio poe' uma
// traducao entre eles.  Traducao e' lugar onde a diferenca medida pode nascer
// -- e uma diferenca que nasce no arquivo de entrada e' indistinguivel de uma
// que nasce no metodo.  Com o mesmo esquema, o MESMO arquivo roda nos dois.
//
// Sao quatro arquivos por caso, como la': `*.par.contr.yaml` (parametros e
// controles), `*.domain.yaml` (geometria), `*.bc.yaml` (contorno) e
// `*.init.yaml` (condicao inicial).  Este leitor junta todos num espaco de
// nomes so', e o caminho de uma chave e' pontilhado:
// `simulation_par.numsteps`.
//
// TRES DECISOES QUE SAO CONTRA O HABITO, e cada uma fecha um modo de falha
// silenciosa que este projeto ja' pagou em outro lugar:
//
//   CHAVE NAO LIDA E' ERRO.  `convecdiscrtype: cubbista`, com um `b` a mais,
//   nao casa com nenhum nome conhecido, cai no padrao e roda O ESQUEMA ERRADO
//   sem dizer nada.  O leitor registra o que foi consultado e sabe dizer, ao
//   final, quais chaves do arquivo ninguem leu;
//
//   NAO HA' PADRAO ESCONDIDO NO LEITOR.  Ou a chave e' obrigatoria -- e faltar
//   lanca, dizendo o caminho completo --, ou quem chama passa o padrao
//   explicitamente na propria chamada, onde quem le' o codigo o ve';
//
//   O TIPO E' COBRADO.  `dt: "0.001"` com aspas e' texto, nao numero, e um
//   leitor que converte por conveniencia aceita `dt: muito pouco` como zero.

#ifndef T8FLOW_COMUM_IO_CONFIGURACAO_HXX
#define T8FLOW_COMUM_IO_CONFIGURACAO_HXX

#include <map>
#include <set>
#include <string>
#include <vector>

namespace t8flow {

class Configuracao {
 public:
  Configuracao() = default;

  //! Le um arquivo e JUNTA ao que ja' existe.
  //
  // Juntar e' aditivo e a colisao e' ERRO, nao sobrescrita: dois arquivos do
  // mesmo caso declarando `simulation_par.dt` com valores diferentes e' defeito
  // do caso, e o ultimo a ser lido vencer esconde isso.
  void junta(const std::string &caminho);

  //! Le os quatro arquivos de um caso pelo prefixo, como o HiGFlow os nomeia:
  //! `<prefixo>.par.contr.yaml`, `.domain.yaml`, `.bc.yaml`, `.init.yaml`.
  //! Os que nao existirem sao PULADOS -- um caso sem condicao inicial e'
  //! legitimo --, e os que existirem e nao abrirem lancam.
  void junta_caso(const std::string &prefixo);

  // ---- leitura obrigatoria -------------------------------------------------
  //
  // Faltar LANCA, e a mensagem traz o caminho pontilhado completo e o arquivo
  // onde o resto do bloco foi encontrado, que e' o que permite achar o erro
  // sem procurar.
  double real(const std::string &chave) const;
  long inteiro(const std::string &chave) const;
  bool booleano(const std::string &chave) const;
  std::string texto(const std::string &chave) const;

  // ---- leitura com padrao EXPLICITO ---------------------------------------
  double real(const std::string &chave, double padrao) const;
  long inteiro(const std::string &chave, long padrao) const;
  bool booleano(const std::string &chave, bool padrao) const;
  std::string texto(const std::string &chave, const std::string &padrao) const;

  //! Um valor de ENUMERACAO: o texto tem de estar na lista, e a mensagem de
  //! erro TRAZ A LISTA.  Devolve o indice.
  //
  // E' aqui que `cubbista` morre: sem a lista, o nome desconhecido viraria o
  // padrao e o esquema trocaria em silencio.
  int opcao(const std::string &chave, const std::vector<std::string> &validos) const;
  int opcao(const std::string &chave, const std::vector<std::string> &validos,
            int padrao) const;

  bool tem(const std::string &chave) const;

  // ---- o que ninguem leu ---------------------------------------------------
  //
  // CHAVE NAO LIDA E' SUSPEITA, e quase sempre e' erro de digitacao ou chave de
  // uma versao antiga do esquema.  Quem monta o caso chama isto depois de
  // configurar tudo e decide o que fazer -- avisar ou recusar.
  std::vector<std::string> nao_lidas() const;

  //! Todas as chaves carregadas, em ordem.  Para diagnostico e para o cabecalho
  //! da saida: relatorio sem a configuracao e' numero sem procedencia.
  std::vector<std::string> chaves() const;

  //! De onde veio cada chave.  Diagnostico de colisao e de caso remendado.
  const std::string &origem(const std::string &chave) const;

 private:
  const std::string *acha(const std::string &chave) const;
  [[noreturn]] void faltou(const std::string &chave) const;

  std::map<std::string, std::string> val_;
  std::map<std::string, std::string> de_;
  mutable std::set<std::string> lidas_;
};

}  // namespace t8flow

#endif  // T8FLOW_COMUM_IO_CONFIGURACAO_HXX
